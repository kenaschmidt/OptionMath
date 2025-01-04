using System.ComponentModel.DataAnnotations;
using System.Diagnostics.Contracts;
using System.Runtime.InteropServices;
using TradingCalendar;

namespace OptionMath
{

    public enum OptionType
    {
        Put = -1,
        Call = 1
    }

    public enum OptionTradingHours
    {
        MarketHours = 0,
        TwentyFourHours = 1
    }

    public static class OptionPricing
    {

        public static int DayCountStandard { get; set; } = 365;


        #region Statistical Functions

        // Cumulative distribution function for the standard normal distribution
        internal static double N(double x)
        {
            double a1 = 0.31938153;
            double a2 = -0.356563782;
            double a3 = 1.781477937;
            double a4 = -1.821255978;
            double a5 = 1.330274429;
            double k = 1 / (1 + 0.2316419 * Math.Abs(x));
            double y = 1 - 1 / Math.Sqrt(2 * Math.PI) * Math.Exp(-Math.Pow(x, 2) / 2)
                       * (a1 * k + a2 * Math.Pow(k, 2) + a3 * Math.Pow(k, 3)
                          + a4 * Math.Pow(k, 4) + a5 * Math.Pow(k, 5));
            return x < 0 ? 1 - y : y;
        }

        // Standardized Normal Density Function n(x)
        internal static double n(double x)
        {
            // As defined on page 353

            var a = 1 / (Math.Sqrt(2 * Math.PI));
            var b = Math.Exp(-Math.Pow(x, 2) / 2);

            return a * b;
        }

        internal static double d1(double underlyingPrice, double optionStrikePrice, double timeToExpiration, double riskFreeRatePercent, double volatilityPercent, double dividendYieldPercent)
        {
            // As defined on page 295

            double So = underlyingPrice;
            double T = timeToExpiration;
            double v = volatilityPercent;
            double r = riskFreeRatePercent;
            double q = dividendYieldPercent;
            double K = optionStrikePrice;

            double a1 = Math.Log(So / K);
            double a2 = (r - q + (Math.Pow(v, 2) / 2)) * T;
            double b = v * Math.Sqrt(T);

            var ret = (a1 + a2) / b;

            return ret;
        }

        internal static double d2(double underlyingPrice, double optionStrikePrice, double timeToExpiration, double riskFreeRatePercent, double volatilityPercent, double dividendYieldPercent)
        {
            // As defined on page 295

            double T = timeToExpiration;
            double v = volatilityPercent;

            var d_2 = d1(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent, dividendYieldPercent) - (v * Math.Sqrt(T));

            return d_2;
        }

        #endregion

        #region Option Pricing Calculations

        /// <summary>
        /// Calculates the price of a European-style option (no early exercise) 
        /// </summary>
        /// <param name="optionType">Put or Call</param>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <returns></returns>
        public static double EuropeanOptionPrice(OptionType optionType, double underlyingPrice, double optionStrikePrice, double timeToExpiration, double riskFreeRatePercent, double volatilityPercent)
        {
            if (optionType == OptionType.Call)
                return _CalculateEuroCallOptionPrice(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent);
            else
                return _CalculateEuroPutOptionPrice(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent);
        }

        /// <summary>
        /// Calculates the price of an American-style option (early exercise allowed)
        /// </summary>
        /// <param name="optionType">Put or Call</param>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="riskFreeRatePercent">Risk free interest rate percentage expressed as a decimal</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <returns></returns>
        public static double AmericanOptionPrice(OptionType optionType, double underlyingPrice, double optionStrikePrice, double timeToExpiration, double riskFreeRatePercent, double volatilityPercent)
        {
            if (optionType == OptionType.Call)
                return _CalculateAmericanCallOptionPrice(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent);
            else
                return _CalculateAmericanPutOptionPrice(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent);
        }

        internal static double _CalculateEuroCallOptionPrice(double underlyingPrice, double strikePrice, double timeToExpiration, double riskFreeRate, double volatility)
        {
            double d1 = (Math.Log(underlyingPrice / strikePrice)
                        + (riskFreeRate + Math.Pow(volatility, 2) / 2) * timeToExpiration)
                        / (volatility * Math.Sqrt(timeToExpiration));

            double d2 = d1 - volatility * Math.Sqrt(timeToExpiration);

            double callPrice = underlyingPrice * N(d1)
                             - strikePrice * Math.Exp(-riskFreeRate * timeToExpiration) * N(d2);

            return callPrice;
        }
        internal static double _CalculateEuroPutOptionPrice(double underlyingPrice, double strikePrice, double timeToExpiration, double riskFreeRate, double volatility)
        {
            double d1 = (Math.Log(underlyingPrice / strikePrice) + (riskFreeRate + volatility * volatility / 2.0) * timeToExpiration) / (volatility * Math.Sqrt(timeToExpiration));
            double d2 = d1 - volatility * Math.Sqrt(timeToExpiration);

            double putOptionPrice = strikePrice * Math.Exp(-riskFreeRate * timeToExpiration) * N(-d2)
                - underlyingPrice * N(-d1);

            return putOptionPrice;
        }
        internal static double _CalculateAmericanCallOptionPrice(double underlyingPrice, double strikePrice, double timeToExpiration, double riskFreeRate, double volatility, int numberOfSteps = 100)
        {
            double deltaT = timeToExpiration / numberOfSteps;
            double up = Math.Exp(volatility * Math.Sqrt(deltaT));
            double down = 1.0 / up;
            double pUp = (Math.Exp(riskFreeRate * deltaT) - down) / (up - down);
            double pDown = 1.0 - pUp;

            double[] underlyingPrices = new double[numberOfSteps + 1];
            double[] optionValues = new double[numberOfSteps + 1];

            for (int i = 0; i <= numberOfSteps; i++)
            {
                underlyingPrices[i] = underlyingPrice * Math.Pow(up, numberOfSteps - i) * Math.Pow(down, i);
                optionValues[i] = Math.Max(underlyingPrices[i] - strikePrice, 0);
            }

            for (int i = numberOfSteps - 1; i >= 0; i--)
            {
                for (int j = 0; j <= i; j++)
                {
                    underlyingPrices[j] = underlyingPrices[j] / up;
                    optionValues[j] = (pUp * optionValues[j] + pDown * optionValues[j + 1]) * Math.Exp(-riskFreeRate * deltaT);
                    optionValues[j] = Math.Max(optionValues[j], underlyingPrices[j] - strikePrice);
                }
            }

            return optionValues[0];
        }
        internal static double _CalculateAmericanPutOptionPrice(double underlyingPrice, double strikePrice, double timeToExpiration, double riskFreeRate, double volatility, int numberOfSteps = 100)
        {
            double deltaT = timeToExpiration / numberOfSteps;
            double up = Math.Exp(volatility * Math.Sqrt(deltaT));
            double down = 1.0 / up;
            double pUp = (Math.Exp(riskFreeRate * deltaT) - down) / (up - down);
            double pDown = 1.0 - pUp;

            double[] underlyingPrices = new double[numberOfSteps + 1];
            double[] optionValues = new double[numberOfSteps + 1];

            for (int i = 0; i <= numberOfSteps; i++)
            {
                underlyingPrices[i] = underlyingPrice * Math.Pow(up, numberOfSteps - i) * Math.Pow(down, i);
                optionValues[i] = Math.Max(strikePrice - underlyingPrices[i], 0);
            }

            for (int i = numberOfSteps - 1; i >= 0; i--)
            {
                for (int j = 0; j <= i; j++)
                {
                    underlyingPrices[j] = underlyingPrices[j] / up;
                    optionValues[j] = (pUp * optionValues[j] + pDown * optionValues[j + 1]) * Math.Exp(-riskFreeRate * deltaT);
                    optionValues[j] = Math.Max(optionValues[j], strikePrice - underlyingPrices[j]);
                }
            }

            return optionValues[0];
        }

        /// <summary>
        /// Calculates the implied volatility of an option from a given price, using iterative deduction
        /// </summary>
        /// <param name="optionType">Put or Call</param>
        /// <param name="underlyingPrice">>Current Spot price of the underlying instrument</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="optionPrice">Price of the option</param>
        /// <returns>Implied volatility percentage of the option expressed as a decimal</returns>
        /// <exception cref="ApplicationException">If the iterative deduction cannot determine IV within the acceptable margin of error</exception>
        public static double ImpliedVolatility(OptionType optionType, double underlyingPrice, double optionStrikePrice, double timeToExpiration, double riskFreeRatePercent, double optionPrice)
        {
            // Initial guess for implied volatility
            double volatilityGuess = 1.00;
            double step = 1.00;

            // Maximum number of iterations
            int maxIterations = 100;

            // Tolerance for convergence
            double tolerance = 0.001;

            // Iteratively solve for the implied volatility

            double guessPrice = AmericanOptionPrice(optionType, underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityGuess);

            int iterations = 0;
            while (Math.Abs(guessPrice - optionPrice) > tolerance && iterations < maxIterations)
            {

                // If the calculated price is too low, raise the IV.  If it is too high, lower.

                if (guessPrice < optionPrice)
                {
                    volatilityGuess += step;
                }
                else
                {
                    step /= 2.0;
                    volatilityGuess -= step;
                }

                guessPrice = AmericanOptionPrice(optionType, underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityGuess);
                iterations++;
            }

            if (iterations == maxIterations)
            {
                throw new ApplicationException("*** ERROR: Implied volatility did not converge ***");
            }

            return volatilityGuess;
        }

        /// <summary>
        /// Calculates the Delta of an option, expressed as a decimal value between 0 and 1 (or -1).
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionType">Put or Call</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <returns>A value between 0.0 and 1.0 (for a call) or 0.0 and -1.0 (for a put)</returns>
        /// <exception cref="ArgumentException"></exception>
        public static double Delta(double underlyingPrice, OptionType optionType, double optionStrikePrice, double timeToExpiration, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent)
        {
            if (optionType == OptionType.Call)
            {
                return N(d1(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent, dividendYieldPercent));
            }
            else if (optionType == OptionType.Put)
            {
                return N(d1(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent, dividendYieldPercent)) - 1;
            }
            else
                throw new ArgumentException();
        }

        /// <summary>
        /// Calculates the Gamma of an option - the change in delta per 1.00 point (dollar) change in the underlying spot price.  Always positive for long options.
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <returns>A value between 0.0 and 1.0</returns>
        public static double Gamma(double underlyingPrice, double optionStrikePrice, double timeToExpiration, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent)
        {
            var T = timeToExpiration;

            var a = n(d1(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent, dividendYieldPercent));
            var b = underlyingPrice * volatilityPercent * Math.Sqrt(T);

            return a / b;
        }

        /// <summary>
        /// Calculates the Theta of an option - the change in option price due to change in time, expressed as a decimal value per 1 day of time.  Always negative for long options.
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionType">Put or Call</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <returns>A negative value representing the expected change in option price per 1 day of time</returns>
        /// <exception cref="ArgumentException"></exception>
        public static double Theta(double underlyingPrice, OptionType optionType, double optionStrikePrice, double timeToExpiration, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent)
        {
            // As defined on page 353

            double So = underlyingPrice;
            double T = timeToExpiration;

            if (optionType == OptionType.Call)
            {
                double a1 = So * (n(d1(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent, dividendYieldPercent))) * volatilityPercent * Math.Exp(-dividendYieldPercent * T);
                double a2 = 2 * Math.Sqrt(T);

                double c1 = (dividendYieldPercent * So * N(d1(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent, dividendYieldPercent)) * Math.Exp(-dividendYieldPercent * T));

                double b1 = riskFreeRatePercent * optionStrikePrice * Math.Exp(-riskFreeRatePercent * T) * N(d2(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent, dividendYieldPercent));

                double ret = -(a1 / a2) + c1 - b1;

                return ret / DayCountStandard;
            }
            else if (optionType == OptionType.Put)
            {
                double a1 = So * n(d1(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent, dividendYieldPercent)) * volatilityPercent * Math.Exp(-dividendYieldPercent * T);
                double a2 = 2 * Math.Sqrt(T);

                double c1 = (dividendYieldPercent * So * N(-d1(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent, dividendYieldPercent)) * Math.Exp(-dividendYieldPercent * T));

                double b1 = riskFreeRatePercent * optionStrikePrice * Math.Exp(-riskFreeRatePercent * T) * N(-d2(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent, dividendYieldPercent));

                double ret = -(a1 / a2) - c1 + b1;

                return ret / DayCountStandard;
            }
            else
                throw new ArgumentException();
        }

        /// <summary>
        /// Calculates the Vega of an option - the change in option price per 1.00% change in option implied volatility.
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <returns>A decimal value representing the expected change in option price per 1.00% change in implied volatility</returns>
        public static double Vega(double underlyingPrice, double optionStrikePrice, double timeToExpiration, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent)
        {
            // As defined on page 361

            double So = underlyingPrice;
            double T = timeToExpiration;

            double d_1 = d1(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent, dividendYieldPercent);

            return (So * Math.Sqrt(T) * n(d_1)) / 100;
        }

        /// <summary>
        /// Calculates the Rho of an option - the change in option price per 1.00% change in the annualized risk free interest rate.
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionType">Put or Call</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <returns></returns>
        /// <exception cref="NotImplementedException"></exception>
        public static double Rho(double underlyingPrice, OptionType optionType, double optionStrikePrice, double timeToExpiration, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent)
        {
            if (optionType == OptionType.Call)
                return ((optionStrikePrice) * (timeToExpiration) * (Math.Exp((-riskFreeRatePercent * timeToExpiration)) * N(d2(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent, dividendYieldPercent))));
            else if (optionType == OptionType.Put)
                return (-(optionStrikePrice) * (timeToExpiration) * (Math.Exp((-riskFreeRatePercent * timeToExpiration)) * N(-d2(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent, dividendYieldPercent))));

            else throw new ArgumentException();
        }

        /// <summary>
        /// Calculates the Lambda of an option - the change in option price per 1.00% change in spot price.  Also known as option elasticity or leverage.
        /// </summary>
        /// <param name="underlyingPrice">Current spot price of the instrument</param>
        /// <param name="optionPrice">Current price of the option</param>
        /// <param name="optionDelta">Current Delta of the option</param>
        /// <returns>Value representing the change in option price</returns>
        public static double Lambda(double underlyingPrice, double optionPrice, double optionDelta)
        {
            // Calculate the percentage change in the underlying price
            double underlyingPercentChange = 0.01;

            // Calculate the change in the option price for a 1% change in the underlying price
            double optionPriceChange = optionDelta * underlyingPrice * underlyingPercentChange;

            // Calculate the Lambda of the option
            double optionLambda = optionPriceChange / optionPrice * 100;

            return optionLambda;
        }

        /// <summary>
        /// Calculates the Vanna of an option - the change in option Delta per 1.00% change in option implied volatility, or the change in option Vega per 1.00 point (dollar) change in the underlying spot price.
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionType">Put or Call</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <returns>Vanna of an option expressed as a decimal</returns>
        public static double Vanna(double underlyingPrice, OptionType optionType, double optionStrikePrice, double timeToExpiration, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent)
        {
            double delta1 = Delta(underlyingPrice, optionType, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
            double delta2 = Delta(underlyingPrice, optionType, optionStrikePrice, timeToExpiration, (volatilityPercent + .01), riskFreeRatePercent, dividendYieldPercent);

            return (delta2 - delta1);
        }

        /// <summary>
        /// Calculates the Charm of an option - the change in option Delta per 1.00 day in time.  Also known as dDelta/dTime (Delta Time Decay)
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionType">Put or Call</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <returns>Charm of an option expressed as a decimal value of delta, normalized to 1.0 day</returns>
        public static double Charm(double underlyingPrice, OptionType optionType, double optionStrikePrice, double timeToExpiration, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent)
        {
            if (timeToExpiration <= (1.0 / DayCountStandard))
            {
                double cDelta = Delta(underlyingPrice, optionType, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
                return -cDelta;
            }

            double delta1 = Delta(underlyingPrice, optionType, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
            double delta2 = Delta(underlyingPrice, optionType, optionStrikePrice, timeToExpiration - (1.0 / DayCountStandard), volatilityPercent, riskFreeRatePercent, dividendYieldPercent);

            return (delta2 - delta1);
        }

        /// <summary>
        /// Calculates the Vomma (or Volga) of an option - the change in option Vega per 1.0% change in option implied volatility
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <returns>A decimal value representing the change in Vega per 1.0% change in implied volatility</returns>
        public static double Vomma(double underlyingPrice, double optionStrikePrice, double timeToExpiration, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent)
        {
            double vega1 = Vega(underlyingPrice, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
            double vega2 = Vega(underlyingPrice, optionStrikePrice, timeToExpiration, (volatilityPercent + .01), riskFreeRatePercent, dividendYieldPercent);

            return (vega2 - vega1);
        }

        /// <summary>
        /// Calculates the Veta for an option - the change in Vega per 1.0 day (Vega Time Decay)
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <returns>The Veta of an option, expressed as a decimal</returns>
        public static double Veta(double underlyingPrice, double optionStrikePrice, double timeToExpiration, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent)
        {
            if (timeToExpiration <= (1.0 / DayCountStandard))
            {
                double cVega = Vega(underlyingPrice, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
                return -cVega;
            }

            double vega1 = Vega(underlyingPrice, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
            double vega2 = Vega(underlyingPrice, optionStrikePrice, (timeToExpiration - (1.0 / DayCountStandard)), volatilityPercent, riskFreeRatePercent, dividendYieldPercent);

            return (vega2 - vega1);
        }

        /// <summary>
        /// Calculates the Speed of an option - the change in option Gamma per 1.0 point (dollar) change in the underlying spot price.
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <returns>Option Speed expressed as a decimal</returns>
        public static double Speed(double underlyingPrice, double optionStrikePrice, double timeToExpiration, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent)
        {
            double gamma = Gamma(underlyingPrice, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
            double d_1 = d1(underlyingPrice, optionStrikePrice, timeToExpiration, riskFreeRatePercent, volatilityPercent, dividendYieldPercent);
            double T = timeToExpiration;

            return -(gamma * (1 + (d_1 / (volatilityPercent * Math.Sqrt(T)))) / underlyingPrice);
        }

        /// <summary>
        /// Calculates the Zomma of an option - the change in option Gamma per 1.0% change in implied volatility
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <returns>Option Zomma expressed as a decimal</returns>
        public static double Zomma(double underlyingPrice, double optionStrikePrice, double timeToExpiration, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent)
        {
            double gamma1 = Gamma(underlyingPrice, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
            double gamma2 = Gamma(underlyingPrice, optionStrikePrice, timeToExpiration, (volatilityPercent + .01), riskFreeRatePercent, dividendYieldPercent);

            return (gamma2 - gamma1);
        }

        /// <summary>
        /// Calculates the Color of an option - the change in option Gamma per 1.0 day of time. (Gamma Time Decay)
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <returns>The change in Gamma per 1.0 day passage of time, expressed as a decimal</returns>
        public static double Color(double underlyingPrice, double optionStrikePrice, double timeToExpiration, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent)
        {
            if (timeToExpiration <= (1.0 / DayCountStandard))
            {
                double cGamma1 = Gamma(underlyingPrice, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
                double cGamma2 = Gamma(underlyingPrice, optionStrikePrice, 0.00001, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);

                return (cGamma2 - cGamma1);
            }

            double gamma1 = Gamma(underlyingPrice, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
            double gamma2 = Gamma(underlyingPrice, optionStrikePrice, (timeToExpiration - (1.0 / DayCountStandard)), volatilityPercent, riskFreeRatePercent, dividendYieldPercent);

            return (gamma2 - gamma1);
        }

        #endregion

        #region Notional Calculations

        //
        // Delta
        //

        /// <summary>
        /// Calculates the notional Delta value of an option position given the number of open contracts and current pricing.
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionType">Put or Call</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="contractCount">Whole number of open contracts, can be positive or negative</param>       
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <param name="contractMultiplier">Number of underlying shares represented by one contract.  Default is 100.</param>
        /// <returns>Net underlying notional dollar value of an option position.</returns>
        /// <exception cref="ArgumentException"></exception>
        public static double Delta_Notional_From_Contracts(double underlyingPrice, OptionType optionType, double optionStrikePrice, double timeToExpiration, int contractCount, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent, double contractMultiplier = 100)
        {
            var deltaDecimal = Delta(underlyingPrice, optionType, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
            return (deltaDecimal * optionStrikePrice * contractMultiplier * contractCount);
        }

        /// <summary>
        /// Calculates the number of open contracts based on a provided notional Delta dollar value.
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionType">Put or Call</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="deltaNotional">Net notional dollar value of an option position</param>       
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <param name="contractMultiplier">Number of underlying shares represented by one contract.  Default is 100.</param>
        /// <returns>Number of open contracts for an option position.</returns>
        /// <exception cref="ArgumentException"></exception>
        public static double Contracts_From_Delta_Notional(double underlyingPrice, OptionType optionType, double optionStrikePrice, double timeToExpiration, double deltaNotional, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent, double contractMultiplier = 100)
        {
            var deltaDecimal = Delta(underlyingPrice, optionType, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
            return (deltaNotional / contractMultiplier / optionStrikePrice / deltaDecimal);
        }

        //
        // Gamma
        //

        /// <summary>
        /// Calculates the notional Gamma value of an option position given the number of open contracts and current pricing.
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="contractCount">Whole number of open contracts, can be positive or negative</param>       
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <param name="contractMultiplier">Number of underlying shares represented by one contract.  Default is 100.</param>
        /// <returns>Net notional Gamma for a position - the change in notional underlying value per 1.00 point (dollar) change in spot price of the underlying.</returns>
        /// <exception cref="ArgumentException"></exception>
        public static double Gamma_Notional_From_Contracts(double underlyingPrice, double optionStrikePrice, double timeToExpiration, double contractCount, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent, double contractMultiplier = 100)
        {
            var gammaDecimal = Gamma(underlyingPrice, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
            return (gammaDecimal * optionStrikePrice * contractMultiplier * contractCount);
        }

        //
        // Theta
        //

        /// <summary>
        /// Calculates the notional Theta value of an option position given the number of open contracts and current pricing.
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionType">Put or Call</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="contractCount">Whole number of open contracts, can be positive or negative</param>       
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <param name="contractMultiplier">Number of underlying shares represented by one contract.  Default is 100.</param>
        /// <returns>Net notional Theta dollar value of an option position - the total dollar amount of decay for a 1.0 day period</returns>
        /// <exception cref="ArgumentException"></exception>
        public static double Theta_Notional_From_Contracts(double underlyingPrice, OptionType optionType, double optionStrikePrice, double timeToExpiration, double contractCount, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent, double contractMultiplier = 100)
        {
            var thetaDecimal = Theta(underlyingPrice, optionType, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
            return (thetaDecimal * contractCount * contractMultiplier);
        }

        //
        // Vega
        //

        /// <summary>
        /// Calculates the notional Vega value of an option position given the number of open contracts and current pricing.
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="contractCount">Whole number of open contracts, can be positive or negative</param>       
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <param name="contractMultiplier">Number of underlying shares represented by one contract.  Default is 100.</param>
        /// <returns>Net notional Vega dollar value of an option position - the total dollar price (value) change given a 1.00% change in implied volatility.</returns>
        /// <exception cref="ArgumentException"></exception>
        public static double Vega_Notional_From_Contracts(double underlyingPrice, double optionStrikePrice, double timeToExpiration, double contractCount, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent, double contractMultiplier = 100)
        {
            var vegaDecimal = Vega(underlyingPrice, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
            return (vegaDecimal * contractCount * contractMultiplier);
        }

        //
        // Rho
        //

        /// <summary>
        /// Calculates the notional Rho value of an option position given the number of open contracts and current pricing.
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>        
        /// <param name="optionType">Put or Call</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="contractCount">Whole number of open contracts, can be positive or negative</param>       
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <param name="contractMultiplier">Number of underlying shares represented by one contract.  Default is 100.</param>
        /// <returns>Net notional Rho dollar value of an option position - the total dollar price (value) change given a 1.00% change in the risk-free interest rate.</returns>
        /// <exception cref="ArgumentException"></exception>
        public static double Rho_Notional_From_Contracts(double underlyingPrice, OptionType optionType, double optionStrikePrice, double timeToExpiration, double contractCount, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent, double contractMultiplier = 100)
        {
            var rhoDecimal = Rho(underlyingPrice, optionType, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
            return (rhoDecimal * contractCount * contractMultiplier);
        }

        //
        // Vanna
        //

        /// <summary>
        /// Calculates the notional Vanna value of an option position given the number of open contracts and current pricing.
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>        
        /// <param name="optionType">Put or Call</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="contractCount">Whole number of open contracts, can be positive or negative</param>       
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <param name="contractMultiplier">Number of underlying shares represented by one contract.  Default is 100.</param>
        /// <returns>Net notional Vanna dollar value of an option position -  the change in notional underlying delta value per 1.00% change in implied volatility.</returns>
        /// <exception cref="ArgumentException"></exception>
        public static double Vanna_Notional_From_Contracts(double underlyingPrice, OptionType optionType, double optionStrikePrice, double timeToExpiration, double contractCount, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent, double contractMultiplier = 100)
        {
            var vannaDecimal = Vanna(underlyingPrice, optionType, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
            return (vannaDecimal * optionStrikePrice * contractMultiplier * contractCount);
        }

        //
        // Charm
        //

        /// <summary>
        /// Calculates the notional Charm value of an option position given the number of open contracts and current pricing.
        /// </summary>
        /// <param name="underlyingPrice">Current Spot price of the underlying instrument</param>        
        /// <param name="optionType">Put or Call</param>
        /// <param name="optionStrikePrice">Strike price of the option</param>
        /// <param name="timeToExpiration">Time to expiration expressed as a fraction of a trading year</param>
        /// <param name="volatilityPercent">Volatility percentage expressed as a decimal</param>
        /// <param name="contractCount">Whole number of open contracts, can be positive or negative</param>       
        /// <param name="riskFreeRatePercent">Risk free interest rate percentagee expressed as a decimal</param>
        /// <param name="dividendYieldPercent">Dividend yield percentage expressed as a decimal</param>
        /// <param name="contractMultiplier">Number of underlying shares represented by one contract.  Default is 100.</param>
        /// <returns>Net notional Charm dollar value of an option position -  the change in notional underlying delta value per 1.0 day of time.</returns>
        /// <exception cref="ArgumentException"></exception>
        public static double Charm_Notional_From_Contracts(double underlyingPrice, OptionType optionType, double optionStrikePrice, double timeToExpiration, double contractCount, double volatilityPercent, double riskFreeRatePercent, double dividendYieldPercent, double contractMultiplier = 100)
        {
            var charmDecimal = Charm(underlyingPrice, optionType, optionStrikePrice, timeToExpiration, volatilityPercent, riskFreeRatePercent, dividendYieldPercent);
            return (charmDecimal * optionStrikePrice * contractMultiplier * contractCount);
        }


        //
        // Vomma
        //

        //
        // Veta
        //

        //
        // Speed
        //

        //
        // Color
        //

        //
        // Zomma
        //

        #endregion
    }
}
