## Copyright 2012 Patrick Steele.
##
## This file is part of curvefit.
##
## Curvefit is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Curvefit is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Curvefit. If not, see <http://www.gnu.org/licenses/>.

library(hash);

curvefit = function(x, ...) {
  UseMethod("curvefit");
};

curvefit.default = function(formula, data, prior_mean, max_knots=NULL, 
  c_param=.4, poly_deg=2, knot_continuity=1,
  mse_relative=10^-3, mse_absolute=10^3, diagnostics=FALSE, ...) {
  #
  # Data should be an S formula of the form y ~ x, where data$y is the
  # response to data$x.
  #
  # Arguments:
  #  max_knots The maximum number of knots in the model.
  #  prior_mean The prior mean of the number of knots.
  #  c_param Controls rate of mixing birth, death, and move steps. Needs
  #    to be in (0, .5).
  #  poly_deg The degree of the polynomials fit between knot points. (l
  #    in the paper)
  #  knot_continuity Controls the degree of continuiting at knot
  #    points. The resulting function will have knotcontinuiting - 1
  #    continuous derivatives at the knot points. (l0 in the paper)
  #

  # Process the formula; coerce xs and ys to be vectors of the data,
  # where ys ~ xs.
  mf = model.frame(formula=formula, data=data);
  xs = model.matrix(attr(mf, "terms"), data=mf);
  ys = model.response(mf);

  data_ndx = 0;
  for (i in 1:length(attr(xs, "assign"))) {
    if (attr(xs, "assign")[i] > 1) {
      stop("Error: Formula contains more than one independent variable;",
          " curvefit can only accept one independent variable.\n");
    }
    else if (attr(xs, "assign")[i] == 1) {
      data_ndx = i;
    }
  }

  xs = as.vector(xs[,data_ndx]);
  ys = as.vector(ys);

  if (length(xs) != length(ys)) {
    stop(paste("length of data vectors must agree. Predictor has length ",
              length(xs),
              ", while response has length ",
              length(ys), sep=""));
  }

  if (is.null(max_knots)) {
    max_knots = 4 * length(xs);
  }

  # Process max_knots.
  max_knots = floor(max_knots);
  if (max_knots < 1) {
    stop(paste("max_knots must be a strictly positive integer, not ",
        max_knots, "\n", sep=""));
  }

  # Check c_param
  if (c_param <= 0 || .5 <= c_param) {
    stop("c_param must be in the interval (0, .5)");
  }

  #
  # We now define a number of useful internal methods.
  #

  # Our prior distribution over the number of knots is a Poisson
  # distribution with mean prior_mean
  p = function(x) {
    return(dpois(x, prior_mean));
  };

  # The likelihood ratio of model m1 to model m2
  L_ratio = function(m1, m2) {
    n = length(xs);

    m1_res = 0;
    m2_res = 0
    for (i in 1:n) {
      m1_res = m1_res + (ys[i] - curvefit.at(m1, xs[i]))^2 / 2;
      m2_res = m2_res + (ys[i] - curvefit.at(m2, xs[i]))^2 / 2;
    }

    m1_res = m1_res * m1$inv_sig_sq / 2;
    m2_res = m2_res * m2$inv_sig_sq / 2;

    m1_arg = n * log(m1$inv_sig_sq) / 2 - m1_res;
    m2_arg = n * log(m2$inv_sig_sq) / 2 - m2_res;
    
    return(exp(m1_arg - m2_arg));
  }

  # Number of candidate knot locations
  Z = function(model) {
    return(2 * (model$l + 1) + model$k * (2 * model$l + 1));
  };

  # Candidate knot locations
  candidates = function(model, rec=FALSE) {
    eligible = c();
    i = 1 + model$l + 1;
    for (knot in model$knots) {
      while (i + model$l < knot) {
        eligible = c(eligible, i);
        i = i + 1;
      }

      i = knot + model$l + 1;
    }
    while (i < max_knots - model$l) {
      eligible = c(eligible, i);
      i = i + 1;
    }

    return(eligible);
  }

  # The domain of the explanatory variable
  dom = c(min(xs), max(xs));

  # Possible knot locations; includes endpoints
  locations = dom[1] + 0:(max_knots + 1) * (dom[2] - dom[1]) / (max_knots + 1);

  # Get a random selection of k0 knots with spacing at least
  # polydeg + 1. We offset by l + 1 since we are sampling over
  # interior knots.
  k0 = rpois(1, prior_mean);
  rand = get_rand_gen(length(locations) - 2 * poly_deg - 2,
    k0, poly_deg + 1);
  selection = rand() + poly_deg + 1;

  # The current model
  model = list();
  model$k = length(selection);
  model$l = poly_deg;
  model$l0 = knot_continuity;
  model$knots = selection;
  model$locations = locations;
  model$formula = get_model_formula(model);

  # Initial coefficient estimate
  data = get_ls_data(model, xs, ys);
  res = lm(model$formula, data);
  model$coef = res$coef;

  # Initial sig_sq
  model$inv_sig_sq = max(10^-10, rgamma(1, shape=10^-3, scale=10^-3));

  if (diagnostics) {
    cat("Interior knots  MSE-absolute  MSE-relative\n");
  }

  #
  # Main loop. We perform birth, death, and move steps until the
  # mean-squared error of the model reaches a steady state. We define
  # 'steady state' to be when the moving average of the change in MSE
  # falls below some absolute level.
  #
  nterms = 10;
  mses = c();
  niters = 0;
  halt = FALSE
  while (halt == FALSE) {

    # Compute birth step, death step, and move step probabilities
    if (model$k == 0) {
      bk = 1;
      dk = 0;
    }
    else {
      bk = c_param * min(1, p(model$k + 1) / p(model$k));
      dk = c_param * min(1, p(model$k - 1) / p(model$k));
    }

    # Form the proposed model
    proposal = model;

    # Make a birth, death, or move step
    u = runif(1);
    birth = FALSE;
    death = FALSE;
    if (u <= bk && Z(model) > 0) {
      # Birth step
      birth = TRUE;
      cs = candidates(proposal);

      if (is.null(cs)) {
        # We cannot grow in this model
        stop(paste("Birth step skipped because max_knots too low; ",
                      "consider increasing", sep=""));
      }
      else {
        added = sample(cs, 1);

        # Update model parameters
        proposal$knots = sort(c(proposal$knots, added));
        proposal$k = proposal$k + 1;
      }
    }
    else if (u <= bk + dk) {
      # Death step
      death = TRUE;
      delete = sample(proposal$knots, 1);
      knots = c();
      for (knot in proposal$knots) {
        if (knot != delete) {
          knots = c(knots, knot);
        }
      }

      # Update model parameters
      proposal$knots = knots;
      proposal$k = proposal$k - 1;
    }
    else {
      # Move step
      cs = candidates(proposal);

      if (is.null(cs)) {
        stop(paste("Move step skipped because max_knots too low; ",
                      "consider increasing", sep=""));
      }
      else {
        move = sample(proposal$knots, 1);
        move_to = sample(cs, 1);

        # Remove the old knot
        knots = c(move_to);
        for (knot in proposal$knots) {
          if (knot != move) {
            knots = c(knots, knot);
          }
        }

        # Update the model parameters
        proposal$knots = sort(knots);
      }
    }

    # Get the polynomial coefficients
    proposal$formula = get_model_formula(proposal);
    data = get_ls_data(proposal, xs, ys);
    res = lm(proposal$formula, data)
    proposal$coef = res$coef;

    # Compute the acceptance probability
    if (birth) {
      prob = L_ratio(proposal, model) * (length(xs) - Z(model)) / length(xs);
      alpha = min(1, prob);
    }
    else if (death) {
      prob = L_ratio(proposal, model) * length(xs) / (length(xs) - Z(model));
      alpha = min(1, prob);
    }
    else {
      alpha = min(1, L_ratio(proposal, model));
    }

    # Accept or reject the proposed model
    if (runif(1) <= alpha) {
      model = proposal;
    }

    #
    # Draw sig_sq using a Gibbs step
    #                                    
    # We know that the full conditional distribution of sig_sq^{-1} is
    # Gamma
    #
    
    sq_residuals = 0;
    for (i in 1:length(xs)) {
      sq_residuals = sq_residuals + (ys[i] - curvefit.at(model, xs[i]))^2;
    }
    sh = 10^-3 + length(xs) / 2;
    sc = (10^-3 + .5 * sq_residuals)^-1;
    model$inv_sig_sq = max(10^-10, rgamma(1, shape=sh, scale=sc));

    # Compute the MSE of the model
    mse = 0;
    for (i in 1:length(xs)) {
      mse = mse + (ys[i] - curvefit.at(model, xs[i]))^2;
    }
    mse = mse / length(xs);

    # Add the mse to the list of the last severals mses
    if (length(mses) == nterms) {
      mses = c(mses[2:length(mses)], mse);
    }
    else {
      mses = c(mses, mse);
    }

    # Check the relative and absolute MSE for the halting condition
    relative = abs(mean(mses) - mse) / mean(mses);
    absolute = mse;
    if (length(mses) >= nterms && relative <= mse_relative
        && absolute <= mse_absolute) {
      halt = TRUE;
    }

    if (diagnostics) {
      cat(paste(model$k, "  ", absolute, "  ", relative, "\n", sep=""));
    }

    # Update the iteration count
    niters = niters + 1;
  }

  # Prepare the final
  model$niters = niters;
  model$mse = mses[length(mses)];

  class(model) = "curvefit";

  return(model);
};

# Evaluates the model using fitted coefficients
curvefit.at = function(model, x) {
  ret = c();

  for (pt in x) {
    res = 0;

    if (pt > model$locations[length(model$locations)]) {
      res = 0
    }
    else {
      for (n in 0 : model$l) {
        r = model$locations[1];
        if (pt >= r) {
          coef = model$coef[[coef_name(n, 0)]];
          if (is.na(coef)) {
            coef = 0;
          }

          if (pt >= r) {
            res = res + coef * (pt - r)^n;
          }
        }
      }

      m = 1;
      while (m <= model$k) {
        r = model$locations[model$knots[m]];
        n = model$l0;
        while (n <= model$l) {
          coef = model$coef[[coef_name(n, m)]];
          if (is.na(coef)) {
            coef = 0;
          }
          if (pt > r) {
            res = res + coef * (pt - r)^n;
          }
          n = n + 1;
        }
        m = m + 1;
      }
    }
    ret = c(ret, res);
  }

  return(ret);
};

print.curvefit = function(x, ...) {
  res = paste(
    "Iterations: ", x$niters, "\n",
    "Knot count: ", x$k + 2, "\n",
    "Mean-squared error of model: ", x$mse, "\n",
    sep="");
  
  cat(res);
};

#
# Utility methods
#

get_rand_gen= function(n, k, l) {
  #
  # Returns a function that will return a random selection of k elements
  # from n elements, where each element is at least l from each other.
  #
  # Example use:
  #
  #  > rand = get_rand_gen(6, 2, 2)
  #  > rand() # Returns c(1, 4), c(1, 5), c(1, 6), c(2, 5), c(2, 6), or c(3, 6)
  #           # with probability 1/6 each
  #

  # Converts an R vector, collection, or numeric to a hash key
  to_key = function(x) {
    if (is.vector(x)) {
      return(paste(as.character(x), collapse=","));
    }
    else {
      return(as.character(x));
    }
  };

  #
  # Count the number of ways to select k elements from n, with spacing
  # l. Use the memo hashtable for memoization.
  #
  memo = hash();
  count = function(n, k) {
    key = to_key(c(n, k));
    if (has.key(key, memo)) {
      return(memo[[key]]);
    }

    # Choosing k elements with l elements between them requires at
    # least (k - 1) * (l + 1) + 1 elements
    if ((k - 1) * (l + 1) + 1 > n) {
      return(0);
    }
    else if ((k - 1) * (l + 1) + 1 == n) {
      return(1);
    }
    else if (k == 0) {
      return(1);
    }
    else if (n < 0) {
      return(0);
    }
    else if (k == 1) {
      return(n);
    }

    else {
      # We can either choose the first element, or not
      choose_first = count(n - 1 - l, k - 1);
      choose_other = count(n - 1, k);
      res = choose_first + choose_other;
      memo[key] = res;
      return(res);
    }
  }
  the_count = count(n, k);
  
  # Pick a random selection

  rand = function() {

    ndx = floor(runif(1) * (the_count + 1));
    if (ndx == the_count + 1) {
      ndx = ndx - 1; # 'zero probability' event that can still happen
                     # when runif(1) == 1; it shouldn't happen, but
                     # the man pages suggest it could under extreme
                     # conditions
    }
    result = c();

    i = 1;
    rem = k;
    while (i <= n && rem > 0) {
      p = count(n - i - l, rem - 1);
      q = count(n - i, rem);
      if (runif(1) * (p + q) <= p) {
        result = c(result, i);
        i = i + 1 + l;
        rem = rem - 1;
      }
      else {
        i = i + 1;
      }
    }

    return(result);
  }
  
  return(rand);
};

#
# This function returns the formula associated with model. y is the
# dependent variable, while the independent variable associated with
# the b_{i, j} coefficient is x_i_j.
#
get_model_formula = function(model) {
  vars = c();
  for (n in 0 : model$l) {
    vars = c(vars, coef_name(n, 0));
  }

  m = 1;
  while (m <= model$k) {
    n = model$l0;
    while (n <= model$l) {
      vars = c(vars, coef_name(n, m));
      n = n + 1;
    }
    m = m + 1;
  }

  return(as.formula(paste("y ~ ", paste(vars, collapse=" + "), " -1", sep="")));
};

#
# Computes the synthetic data for least-squares regression. We
# associate a new variable with each b_{i, j} coefficient to be
# fit. Returns a list of the modified values, indexed by the variable
# name.
#
get_ls_data = function(model, xs, ys) {
  data = list();

  for (n in 0 : model$l) {
    key = coef_name(n, 0);
    cdata = c();

    r = model$locations[1];
    for (x in xs) {
      if (x > r) {
        cdata = c(cdata, (x - r)^n);
      }
      else {
        cdata = c(cdata, 0);
      }
    }

    data[[key]] = cdata;
  }

  m = 1;
  while (m <= model$k) {
    r = model$locations[model$knots[m]];
    n = model$l0;
    while (n <= model$l) {
      key = coef_name(n, m);
      cdata = c();

      for (x in xs) {
        if (x > r) {
          cdata = c(cdata, (x - r)^n);
        }
        else {
          cdata = c(cdata, 0);
        }
      }

      data[[key]] = cdata;
      n = n + 1;
    }
    m = m + 1;
  }

  # Attach response data
  data$y = ys;

  return(data);
}

#
# Returns the name of the independent variable we create for
# coefficient b_{i, j}.
#
coef_name = function(n, m) {
  return(paste("x_", n, "_", m, sep=""));
}
