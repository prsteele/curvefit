source("curvefit.r");

data(package="datasets")

x = 1:200
y = EuStockMarkets[x, 1]

fit = curvefit(y ~ x, list(y=y, x=x), prior_mean=10, mse_relative=10^-2)

plot(x, curvefit.at(fit, x), type="l", color="red")
points(x, y, type="l", color="blue");

## library(datasets)

## nterms = 200;

## xs = 1:nterms;
## ys = EuStockMarkets[xs,1];

## result = curvefit(ys ~ xs, list(xs=xs, ys=ys), poly_deg = 4,
##   knot_continuity=1, prior_mean=20, max_knots=200);

## step = .1;

## yfits = c();
## xfits = c();

## x = min(xs);
## i = 1;
## while (x <= max(xs)) {
##   xfits = c(xfits, x);
##   yfits = c(yfits, curvefit.at(result, x));
##   x = x + step;
## }

## plot(xfits, yfits, type="l", pch=".", col="blue");
## points(xs, ys, type="l", col="red");

## knots_x = 0 * (1 : length(result$knots));
## knots_y = 0 * (1 : length(result$knots));
## for (i in 1:length(result$knots)) {
##   knots_x[i] = result$location[result$knots[i]];
##   knots_y[i] = curvefit.at(result, knots_x[i]);
## }

## points(knots_x, knots_y, type="p", col="blue");

## xs = 1:10;
## ys = c(0 * (1:5), 1 + 0 * (6:10));

## res = curvefit(ys ~ xs, list(ys=ys, xs=xs), poly_deg=1, knot_continuity=0, prior_mean=1, max_knots=10);

## step = .1;

## yfits = c();
## xfits = c();

## x = min(xs);
## i = 1;
## while (x <= max(xs)) {
##   xfits = c(xfits, x);
##   yfits = c(yfits, curvefit.at(res, x));
##   x = x + step;
## }

## plot(xfits, yfits, type="l", pch=".", col="blue");
## points(xs, ys, type="l", col="red");



