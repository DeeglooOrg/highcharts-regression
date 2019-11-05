(function(factory) {
  "use strict";

  if (typeof module === "object" && module.exports) {
    module.exports = factory;
  } else {
    if (typeof define === "function" && define.amd) {
      define(function() {
        return factory;
      });
    } else {
      if (typeof Highcharts !== "undefined") {
        factory(Highcharts);
      } else {
        void 0;
      }
    }
  }
})(function(H) {
  var processSerie = function(s, method, chart) {
    if (s.regression && !s.rendered) {
      s.regressionSettings = s.regressionSettings || {};
      s.regressionSettings.tooltip = s.regressionSettings.tooltip || {};
      s.regressionSettings.dashStyle =
        s.regressionSettings.dashStyle || "solid";
      s.regressionSettings.decimalPlaces =
        s.regressionSettings.decimalPlaces || 2;
      s.regressionSettings.useAllSeries =
        s.regressionSettings.useAllSeries || false;

      var regressionType = s.regressionSettings.type || "polynomial";
      var regression;
      var extraSerie = {
        data: [],
        color: s.regressionSettings.color || (s.color + "66") || "",
        yAxis: s.yAxis,
        lineWidth: s.regressionSettings.lineWidth || 2,
        marker: {
          enabled:
            typeof s.regressionSettings.marker == "undefined"
              ? false
              : s.regressionSettings.marker.enabled,
          symbol:
            typeof s.regressionSettings.marker == "undefined"
              ? undefined
              : s.regressionSettings.marker.symbol,
          radius:
            typeof s.regressionSettings.marker == "undefined"
              ? 4
              : s.regressionSettings.marker.radius
        },
        isRegressionLine: true,
        visible: s.regressionSettings.visible,
        type: s.regressionSettings.linetype || "spline",
        name: s.regressionSettings.name || `${s.name} Trendline`,
        id: s.regressionSettings.id,
        dashStyle: s.regressionSettings.dashStyle || "solid",
        showInLegend: false, //!s.regressionSettings.hideInLegend,
        tooltip: s.regressionSettings.tooltip
      };

      extraSerie.dataLabels = s.regressionSettings.dataLabels || {}
      extraSerie.dataLabels.color = extraSerie.color
      extraSerie.dataLabels.borderColor = extraSerie.borderColor || extraSerie.color

      if (typeof s.regressionSettings.index !== "undefined") {
        extraSerie.index = s.regressionSettings.index;
      }
      if (typeof s.regressionSettings.legendIndex !== "undefined") {
        extraSerie.legendIndex = s.regressionSettings.legendIndex;
      }

      var mergedData = s.data || []
      if (s.regressionSettings.useAllSeries) {
        mergedData = [];
        for (di = 0; di < series.length; di++) {
          var seriesToMerge = series[di];
          mergedData = mergedData.concat(seriesToMerge.data);
        }
      }

      if (regressionType == "linear") {
        var extrapolate = s.regressionSettings.extrapolate || 0;
        regression = _linear(
          mergedData,
          s.regressionSettings.decimalPlaces,
          extrapolate
        );
        extraSerie.type = "line";
      } else if (regressionType == "exponential") {
        var extrapolate = s.regressionSettings.extrapolate || 0;
        regression = _exponential(mergedData, extrapolate);
      } else if (regressionType == "polynomial") {
        const MAX_POLYNOMIAL_ORDER = 6,
          POLYNOMIAL_ORDER_SCALER = 1.2;

        if (mergedData.length == 0) {
          regression = { equation: [], points: [], string: "" };
        } else {
          //Equation for the polynomial order
          var order =
            s.regressionSettings.order ||
            Math.min(
              Math.floor(
                Math.log(POLYNOMIAL_ORDER_SCALER * mergedData.length)
              ) + 1,
              MAX_POLYNOMIAL_ORDER
            );
          var extrapolate = s.regressionSettings.extrapolate || 0;
          regression = _polynomial(mergedData, order, extrapolate);
        }
      } else if (regressionType == "power") {
        var extrapolate = s.regressionSettings.extrapolate || 0;
        regression = _power(mergedData, extrapolate);
      } else if (regressionType == "logarithmic") {
        var extrapolate = s.regressionSettings.extrapolate || 0;
        regression = _logarithmic(mergedData, extrapolate);
      } else if (regressionType == "loess") {
        var loessSmooth = s.regressionSettings.loessSmooth || 25;
        regression = _loess(mergedData, loessSmooth / 100);
      } else {
        console.error("Invalid regression type: ", regressionType);
        return;
      }

      regression.rSquared = coefficientOfDetermination(
        mergedData,
        regression.points
      );
      regression.rValue = _round(
        Math.sqrt(regression.rSquared),
        s.regressionSettings.decimalPlaces
      );
      regression.rSquared = _round(
        regression.rSquared,
        s.regressionSettings.decimalPlaces
      );
      regression.standardError = _round(
        standardError(mergedData, regression.points),
        s.regressionSettings.decimalPlaces
      );
      extraSerie.data = regression.points;
      extraSerie.name = extraSerie.name.replace("%r2", regression.rSquared);
      extraSerie.name = extraSerie.name.replace("%r", regression.rValue);
      extraSerie.name = extraSerie.name.replace(
        "%se",
        regression.standardError
      );

      if (extraSerie.visible === false || regression.string.includes("NaN")) {
        extraSerie.visible = false;
      }
      extraSerie.regressionOutputs = regression;
      return extraSerie;
    }
  };

  H.wrap(H.Chart.prototype, 'update', function (proceed) {
    // console.log("============ Highcharts update ============");
    const regressionSeries = this.series.filter(s => s.options.isRegressionLine);
    const series = this.series.filter(s => !s.options.isRegressionLine);
    const extraSeries = [];
    let i = 0;
    if (series) {
      for (i = 0; i < series.length; i++) {
        const s = series[i];
        if (s.options.regression) {
          const serOptions = (arguments[1].series || []).find(x => x.id === s.options.id) || s.options
          var extraSerie = processSerie(serOptions, 'init', this);
          extraSerie.mdd = JSON.parse(JSON.stringify(serOptions.mdd))
          extraSerie.mdd.label = serOptions.mdd.label && serOptions.mdd.label + " Trendline"
          extraSerie.mddFromData = serOptions.mddFromData && JSON.parse(JSON.stringify(serOptions.mddFromData))
          extraSeries.push(extraSerie);
        }
      }
    }
    if (extraSeries.length) {
      extraSeries.forEach(regressionSerie => {
        const realRegSer = regressionSeries.find(rs => rs.options.id === regressionSerie.id)
        if (realRegSer) {
          realRegSer.update(regressionSerie)
        } else {
          this.addSeries(regressionSerie)
        }
      })
      arguments[1].series = (arguments[1].series || []).concat(extraSeries)
    }
    proceed.apply(this, Array.prototype.slice.call(arguments, 1));
  });

  H.wrap(H.Series.prototype, "setData", function(proceed) {
    // console.log("============ Highcharts setData series ============");
    if (this.options.isRegressionLine) {
      const regLineId = this.options.id
      const regSerie = this.chart.series.find(serie => serie.options.regressionSettings && serie.options.regressionSettings.id === regLineId)
      const extraSerie = processSerie(regSerie.options, 'none', this);
      const newData = extraSerie.data.map((e, i) => {
        return [e[0], e[1], regSerie.options.data && regSerie.options.data[i].date];
      });
      arguments[1] = newData;
    }
    return proceed.apply(this, Array.prototype.slice.call(arguments, 1));
  });

  H.wrap(H.Series.prototype, "drawGraph", function(proceed) {
    // console.log("============ Highcharts drawGraph series ============");
    if (this.options.isRegressionLine) {
      const dates = this.options.data.map(d => d[2])
      this.data = this.data.map((d,i) => {
        d.date = dates[i];
        return d;
      })
    }
    return proceed.apply(this, Array.prototype.slice.call(arguments, 1));
  });
  
  /**
   * Code extracted from https://github.com/Tom-Alexander/regression-js/
   */
  function _exponential(data, extrapolate) {
    var sum = [0, 0, 0, 0, 0, 0],
      n = 0,
      results = [];

    for (len = data.length; n < len; n++) {
      if (data[n].x != null && data[n].y != null && data[n].y > 0) {
        sum[0] += data[n].x; // X
        sum[1] += data[n].y; // Y
        sum[2] += data[n].x * data[n].x * data[n].y; // XXY
        sum[3] += data[n].y * Math.log(data[n].y); // Y Log Y
        sum[4] += data[n].x * data[n].y * Math.log(data[n].y); //YY Log Y
        sum[5] += data[n].x * data[n].y; //XY
      }
    }

    var denominator = sum[1] * sum[2] - sum[5] * sum[5];
    var A = Math.pow(Math.E, (sum[2] * sum[3] - sum[5] * sum[4]) / denominator);
    var B = (sum[1] * sum[4] - sum[5] * sum[3]) / denominator;

    var resultLength = data.length + extrapolate;
    var step = data[data.length - 1][0] - data[data.length - 2][0];

    for (var i = 0, len = resultLength; i < len; i++) {
      var answer = 0;
      if (typeof data[i] !== "undefined") {
        var x = data[i][0];
      } else {
        var x = data[data.length - 1][0] + (i - data.length) * step;
      }

      var coordinate = [x, A * Math.pow(Math.E, B * x)];
      results.push(coordinate);
    }

    results.sort(function(a, b) {
      if (a[0] > b[0]) {
        return 1;
      }
      if (a[0] < b[0]) {
        return -1;
      }
      return 0;
    });

    var string =
      "y = " +
      Math.round(A * 100) / 100 +
      "e^(" +
      Math.round(B * 100) / 100 +
      "x)";

    return { equation: [A, B], points: results, string: string };
  }

  /**
   * Code extracted from https://github.com/Tom-Alexander/regression-js/
   * Human readable formulas:
   *
   *              N * Σ(XY) - Σ(X)
   * intercept = ---------------------
   *              N * Σ(X^2) - Σ(X)^2
   *
   * correlation = N * Σ(XY) - Σ(X) * Σ (Y) / √ (  N * Σ(X^2) - Σ(X) ) * ( N * Σ(Y^2) - Σ(Y)^2 ) ) )
   *
   */
  function _linear(data, decimalPlaces, extrapolate) {
    var sum = [0, 0, 0, 0, 0],
      n = 0,
      results = [],
      N = data.length;

    for (; n < data.length; n++) {
      if (data[n]["x"] != null && data[n].y != null) {
        sum[0] += data[n].x; //Σ(X)
        sum[1] += data[n].y; //Σ(Y)
        sum[2] += data[n].x * data[n].x; //Σ(X^2)
        sum[3] += data[n].x * data[n].y; //Σ(XY)
        sum[4] += data[n].y * data[n].y; //Σ(Y^2)
      } else {
        N -= 1;
      }
    }

    var gradient =
      (N * sum[3] - sum[0] * sum[1]) / (N * sum[2] - sum[0] * sum[0]);
    var intercept = sum[1] / N - (gradient * sum[0]) / N;
    // var correlation = (N * sum[3] - sum[0] * sum[1]) / Math.sqrt((N * sum[2] - sum[0] * sum[0]) * (N * sum[4] - sum[1] * sum[1]));

    var resultLength = data.length + extrapolate;
    var step = data[data.length - 1][0] - data[data.length - 2][0];

    for (var i = 0, len = resultLength; i < len; i++) {
      var answer = 0;
      if (typeof data[i] !== "undefined") {
        var x = data[i][0];
      } else {
        var x = data[data.length - 1][0] + (i - data.length) * step;
      }

      var coorY = x * gradient + intercept;
      if (decimalPlaces) coorY = parseFloat(coorY.toFixed(decimalPlaces));
      var coordinate = [x, coorY];
      results.push(coordinate);
    }

    results.sort(function(a, b) {
      if (a.x > b.x) {
        return 1;
      }
      if (a.x < b.x) {
        return -1;
      }
      return 0;
    });

    var string =
      "y = " +
      Math.round(gradient * 100) / 100 +
      "x + " +
      Math.round(intercept * 100) / 100;
    return { equation: [gradient, intercept], points: results, string: string };
  }

  /**
   *  Code extracted from https://github.com/Tom-Alexander/regression-js/
   */
  function _logarithmic(data, extrapolate) {
    var sum = [0, 0, 0, 0],
      n = 0,
      results = [],
      mean = 0;

    for (len = data.length; n < len; n++) {
      if (data[n].x != null && data[n].y != null && data[n].x > 0) {
        sum[0] += Math.log(data[n].x);
        sum[1] += data[n].y * Math.log(data[n].x);
        sum[2] += data[n].y;
        sum[3] += Math.pow(Math.log(data[n].x), 2);
      }
    }

    var B = (n * sum[1] - sum[2] * sum[0]) / (n * sum[3] - sum[0] * sum[0]);
    var A = (sum[2] - B * sum[0]) / n;

    var resultLength = data.length + extrapolate;
    var step = data[data.length - 1][0] - data[data.length - 2][0];

    for (var i = 0, len = resultLength; i < len; i++) {
      if (typeof data[i] !== "undefined") {
        if (data[i][0] <= 0) {
          continue;
        } else {
          var x = data[i][0];
        }
      } else {
        var x = data[data.length - 1][0] + (i - data.length) * step;
      }

      var coordinate = [x, A + B * Math.log(x)];
      results.push(coordinate);
    }

    results.sort(function(a, b) {
      if (a.x > b.x) {
        return 1;
      }
      if (a.x < b.x) {
        return -1;
      }
      return 0;
    });

    var string =
      "y = " +
      Math.round(A * 100) / 100 +
      " + " +
      Math.round(B * 100) / 100 +
      " ln(x)";

    return { equation: [A, B], points: results, string: string };
  }

  /**
   * Code extracted from https://github.com/Tom-Alexander/regression-js/
   */
  function _power(data, extrapolate) {
    var sum = [0, 0, 0, 0],
      n = 0,
      results = [];

    for (len = data.length; n < len; n++) {
      if (data[n].x != null && data[n].y != null && data[n].x > 0) {
        sum[0] += Math.log(data[n].x);
        sum[1] += Math.log(data[n].y) * Math.log(data[n].x);
        sum[2] += Math.log(data[n].y);
        sum[3] += Math.pow(Math.log(data[n].x), 2);
      }
    }

    var B = (n * sum[1] - sum[2] * sum[0]) / (n * sum[3] - sum[0] * sum[0]);
    var A = Math.pow(Math.E, (sum[2] - B * sum[0]) / n);

    var resultLength = data.length + extrapolate;
    var step = data[data.length - 1][0] - data[data.length - 2][0];

    for (var i = 0, len = resultLength; i < len; i++) {
      var answer = 0;
      if (typeof data[i] !== "undefined") {
        var x = data[i][0];
      } else {
        var x = data[data.length - 1][0] + (i - data.length) * step;
      }

      var coordinate = [x, A * Math.pow(x, B)];
      results.push(coordinate);
    }

    results.sort(function(a, b) {
      if (a.x > b.x) {
        return 1;
      }
      if (a.x < b.x) {
        return -1;
      }
      return 0;
    });

    var string =
      "y = " + Math.round(A * 100) / 100 + "x^" + Math.round(B * 100) / 100;

    return { equation: [A, B], points: results, string: string };
  }

  /**
   * Code extracted from https://github.com/Tom-Alexander/regression-js/
   */
  function _polynomial(data, order, extrapolate) {
    const copiedData = JSON.parse(JSON.stringify(data));
    let firstNonNaNItem = copiedData.find(item => !!item.y);
    let lastYValue = firstNonNaNItem && firstNonNaNItem.y || 0;
    copiedData.forEach(item => !item.y ? (item.y = lastYValue) : (lastYValue = item.y))
    if (typeof order == "undefined") {
      order = 2;
    } else if (order < 0) {
      return {
        equation: [],
        points: copiedData,
        string: "Error! Polynomial order smaller than 0!"
      };
    }

    if (copiedData.length == 1) {
      return {
        equation: [copiedData[0].y],
        points: copiedData,
        string: copiedData[0].y.toString()
      };
    }

    var xAxisCopy = new Array();

    //Normalize input data
    const SCALING_RANGE = 10;
    var max = copiedData.reduce(
      (max, point) => (point.x > max ? point.x : max),
      copiedData[0].x
    );
    var min = copiedData.reduce(
      (min, point) => (point.x < min ? point.x : min),
      copiedData[0].x
    );
    var ratio = SCALING_RANGE / (max - min);

    for (var i = 0; i < copiedData.length; i++) {
      xAxisCopy.push(copiedData[i].x);
      copiedData[i].x = (copiedData[i].x - min) * ratio;
    }

    var lhs = [],
      rhs = [],
      results = [],
      a = 0,
      b = 0,
      i = 0,
      k = order + 1;

    for (; i < k; i++) {
      for (var l = 0, len = copiedData.length; l < len; l++) {
        if (copiedData[l].x != null && copiedData[l].y != null) {
          a += Math.pow(copiedData[l].x, i) * copiedData[l].y;
        }
      }

      lhs.push(a);
      a = 0;
      var c = [];
      for (var j = 0; j < k; j++) {
        for (var l = 0, len = copiedData.length; l < len; l++) {
          if (copiedData[l].y) {
            b += Math.pow(copiedData[l].x, i + j);
          }
        }
        c.push(b);
        b = 0;
      }
      rhs.push(c);
    }
    rhs.push(lhs);

    var equation = gaussianElimination(rhs, k);

    var resultLength = copiedData.length + extrapolate;
    var step = copiedData[copiedData.length - 1].x - copiedData[copiedData.length - 2].x;
    for (var i = 0, len = resultLength; i < len; i++) {
      var answer = 0;
      var x = 0;
      if (typeof copiedData[i] !== "undefined") {
        x = copiedData[i].x;
      } else {
        x = copiedData[copiedData.length - 1].x + (i - copiedData.length) * step;
      }

      for (var w = 0; w < equation.length; w++) {
        answer += equation[w] * Math.pow(x, w);
      }
      results.push([x, answer]);
    }

    var string = "y = ";

    for (var i = equation.length - 1; i >= 0; i--) {
      if (i > 1)
        string += Math.round(equation[i] * 100) / 100 + "x^" + i + " + ";
      else if (i == 1)
        string += Math.round(equation[i] * 100) / 100 + "x" + " + ";
      else string += Math.round(equation[i] * 100) / 100;
    }

    //return the old x-axis values back
    for (var i = 0; i < results.length; i++) {
      results[i][0] = xAxisCopy[i];
      copiedData[i].x = xAxisCopy[i];
    }

    results.sort(function(a, b) {
      if (a.x > b.x) {
        return 1;
      }
      if (a.x < b.x) {
        return -1;
      }
      return 0;
    });

    //WARNING! Results are correct but now the final equation is wrong because it is calculated on the normalized data and then translated back.
    //True equation can be calculated but it is unnecessary cumputation in our case and in most cases you would get really, really small coefficients.

    return { equation: equation, points: results, string: string };
  }

  /**
   * @author: Ignacio Vazquez
   * Based on
   * - http://commons.apache.org/proper/commons-math/download_math.cgi LoesInterpolator.java
   * - https://gist.github.com/avibryant/1151823
   */
  function _loess(data, bandwidth) {
    bandwidth = bandwidth || 0.25;

    var xval = data.map(function(pair) {
      return pair.x || pair[0];
    });
    var distinctX = array_unique(xval);
    if (2 / distinctX.length > bandwidth) {
      bandwidth = Math.min(2 / distinctX.length, 1);
      console.warn("updated bandwith to " + bandwidth);
    }

    var yval = data.map(function(pair) {
      return pair.y || pair[1];
    });

    function array_unique(values) {
      var o = {},
        i,
        l = values.length,
        r = [];
      for (i = 0; i < l; i += 1) o[values[i]] = values[i];
      for (i in o) r.push(o[i]);
      return r;
    }

    function tricube(x) {
      var tmp = 1 - x * x * x;
      return tmp * tmp * tmp;
    }

    var res = [];

    var left = 0;
    var right = Math.floor(bandwidth * xval.length) - 1;

    for (var i in xval) {
      var x = xval[i];

      if (i > 0) {
        if (
          right < xval.length - 1 &&
          xval[right + 1] - xval[i] < xval[i] - xval[left]
        ) {
          left++;
          right++;
        }
      }
      //console.debug("left: "+left  + " right: " + right );
      var edge;
      if (xval[i] - xval[left] > xval[right] - xval[i]) edge = left;
      else edge = right;
      var denom = Math.abs(1.0 / (xval[edge] - x));
      var sumWeights = 0;
      var sumX = 0,
        sumXSquared = 0,
        sumY = 0,
        sumXY = 0;

      var k = left;
      while (k <= right) {
        var xk = xval[k];
        var yk = yval[k];
        var dist;
        if (k < i) {
          dist = x - xk;
        } else {
          dist = xk - x;
        }
        var w = tricube(dist * denom);
        var xkw = xk * w;
        sumWeights += w;
        sumX += xkw;
        sumXSquared += xk * xkw;
        sumY += yk * w;
        sumXY += yk * xkw;
        k++;
      }

      var meanX = sumX / sumWeights;
      //console.debug(meanX);
      var meanY = sumY / sumWeights;
      var meanXY = sumXY / sumWeights;
      var meanXSquared = sumXSquared / sumWeights;

      var beta;
      if (meanXSquared == meanX * meanX) beta = 0;
      else beta = (meanXY - meanX * meanY) / (meanXSquared - meanX * meanX);

      var alpha = meanY - beta * meanX;
      res[i] = beta * x + alpha;
    }
    //console.debug(res);
    return {
      equation: "",
      points: xval.map(function(x, i) {
        return [x, res[i]];
      }),
      string: ""
    };
  }

  /**
   * Code extracted from https://github.com/Tom-Alexander/regression-js/
   */
  function gaussianElimination(a, o) {
    var i = 0,
      j = 0,
      k = 0,
      maxrow = 0,
      tmp = 0,
      n = a.length - 1,
      x = new Array(o);
    for (i = 0; i < n; i++) {
      maxrow = i;
      for (j = i + 1; j < n; j++) {
        if (Math.abs(a[i][j]) > Math.abs(a[i][maxrow])) maxrow = j;
      }
      for (k = i; k < n + 1; k++) {
        tmp = a[k][i];
        a[k][i] = a[k][maxrow];
        a[k][maxrow] = tmp;
      }
      for (j = i + 1; j < n; j++) {
        for (k = n; k >= i; k--) {
          a[k][j] -= (a[k][i] * a[i][j]) / a[i][i];
        }
      }
    }
    for (j = n - 1; j >= 0; j--) {
      tmp = 0;
      for (k = j + 1; k < n; k++) tmp += a[k][j] * x[k];
      x[j] = (a[n][j] - tmp) / a[j][j];
    }
    return x;
  }

  /**
   * @author Ignacio Vazquez
   * See http://en.wikipedia.org/wiki/Coefficient_of_determination for theaorical details
   */
  function coefficientOfDetermination(data, pred) {
    // Sort the initial data { pred array (model's predictions) is sorted  }
    // The initial data must be sorted in the same way in order to calculate the coefficients
    data.sort(function(a, b) {
      if (a[0] > b[0]) {
        return 1;
      }
      if (a[0] < b[0]) {
        return -1;
      }
      return 0;
    });

    // Calc the mean
    var mean = 0;
    var N = data.length;
    for (var i = 0; i < data.length; i++) {
      if (data[i][1] != null) {
        mean += data[i][1];
      } else {
        N--;
      }
    }
    mean /= N;

    // Calc the coefficent of determination
    var SSE = 0;
    var SSYY = 0;
    for (var i = 0; i < data.length; i++) {
      if (data[i][1] != null) {
        SSYY += Math.pow(data[i][1] - pred[i][1], 2);
        SSE += Math.pow(data[i][1] - mean, 2);
      }
    }
    return 1 - SSYY / SSE;
  }

  function standardError(data, pred) {
    var SE = 0,
      N = data.length;

    for (var i = 0; i < data.length; i++) {
      if (data[i][1] != null) {
        SE += Math.pow(data[i][1] - pred[i][1], 2);
      } else {
        N--;
      }
    }
    SE = Math.sqrt(SE / (N - 2));

    return SE;
  }

  function _round(number, decimalPlaces) {
    var decimalFactor = Math.pow(10, decimalPlaces);
    return Math.round(number * decimalFactor) / decimalFactor;
  }
});
