using System;
using System.Diagnostics;
using OptionMath;

namespace OptionMath
{
    public class PerformanceBenchmark
    {
        public static void RunBenchmarks()
        {
            Console.WriteLine("=== OptionMath Performance Benchmark ===\n");

            // Test parameters
            double underlyingPrice = 100.0;
            double strikePrice = 100.0;
            double timeToExpiration = 0.25; // 3 months
            double volatility = 0.25; // 25%
            double riskFreeRate = 0.05; // 5%
            double dividendYield = 0.02; // 2%
            OptionType optionType = OptionType.Call;

            int iterations = 100000;

            // Benchmark 1: Individual Greeks vs Batch Calculation
            Console.WriteLine("Benchmark 1: Individual Greeks vs Batch Calculation");
            Console.WriteLine(new string('-', 60));

            var sw = Stopwatch.StartNew();
            for (int i = 0; i < iterations; i++)
            {
                var delta = OptionPricing.Delta(underlyingPrice, optionType, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
                var gamma = OptionPricing.Gamma(underlyingPrice, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
                var theta = OptionPricing.Theta(underlyingPrice, optionType, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
                var vega = OptionPricing.Vega(underlyingPrice, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
                var rho = OptionPricing.Rho(underlyingPrice, optionType, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
            }
            sw.Stop();
            double individualTime = sw.Elapsed.TotalMilliseconds;
            Console.WriteLine($"Individual Greeks (5 calls): {individualTime:F2} ms");

            sw.Restart();
            for (int i = 0; i < iterations; i++)
            {
                var greeks = OptionPricing.CalculateFirstOrderGreeks(underlyingPrice, optionType, strikePrice,
                    timeToExpiration, volatility, riskFreeRate, dividendYield);
            }
            sw.Stop();
            double batchTime = sw.Elapsed.TotalMilliseconds;
            Console.WriteLine($"Batch Calculation:           {batchTime:F2} ms");
            Console.WriteLine($"Speedup:                     {individualTime / batchTime:F2}x faster");
            Console.WriteLine($"Improvement:                 {(1 - batchTime / individualTime) * 100:F1}%\n");

            // Benchmark 2: All Greeks
            Console.WriteLine("Benchmark 2: All Greeks (13 total)");
            Console.WriteLine(new string('-', 60));

            sw.Restart();
            for (int i = 0; i < iterations / 10; i++)
            {
                var price = OptionPricing.EuropeanOptionPrice(optionType, underlyingPrice, strikePrice, timeToExpiration, riskFreeRate, volatility);
                var delta = OptionPricing.Delta(underlyingPrice, optionType, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
                var gamma = OptionPricing.Gamma(underlyingPrice, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
                var theta = OptionPricing.Theta(underlyingPrice, optionType, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
                var vega = OptionPricing.Vega(underlyingPrice, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
                var rho = OptionPricing.Rho(underlyingPrice, optionType, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
                var vanna = OptionPricing.Vanna(underlyingPrice, optionType, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
                var charm = OptionPricing.Charm(underlyingPrice, optionType, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
                var vomma = OptionPricing.Vomma(underlyingPrice, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
                var veta = OptionPricing.Veta(underlyingPrice, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
                var speed = OptionPricing.Speed(underlyingPrice, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
                var zomma = OptionPricing.Zomma(underlyingPrice, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
                var color = OptionPricing.Color(underlyingPrice, strikePrice, timeToExpiration, volatility, riskFreeRate, dividendYield);
            }
            sw.Stop();
            double allIndividualTime = sw.Elapsed.TotalMilliseconds;
            Console.WriteLine($"Individual calls (13):       {allIndividualTime:F2} ms");

            sw.Restart();
            for (int i = 0; i < iterations / 10; i++)
            {
                var allGreeks = OptionPricing.CalculateAllGreeks(underlyingPrice, optionType, strikePrice,
                    timeToExpiration, volatility, riskFreeRate, dividendYield);
            }
            sw.Stop();
            double allBatchTime = sw.Elapsed.TotalMilliseconds;
            Console.WriteLine($"Batch calculation:           {allBatchTime:F2} ms");
            Console.WriteLine($"Speedup:                     {allIndividualTime / allBatchTime:F2}x faster");
            Console.WriteLine($"Improvement:                 {(1 - allBatchTime / allIndividualTime) * 100:F1}%\n");

            // Benchmark 3: American Options (Binomial Tree)
            Console.WriteLine("Benchmark 3: American Option Pricing (Binomial Tree)");
            Console.WriteLine(new string('-', 60));

            sw.Restart();
            for (int i = 0; i < 1000; i++)
            {
                var americanCall = OptionPricing.AmericanOptionPrice(OptionType.Call, underlyingPrice, strikePrice,
                    timeToExpiration, riskFreeRate, volatility);
                var americanPut = OptionPricing.AmericanOptionPrice(OptionType.Put, underlyingPrice, strikePrice,
                    timeToExpiration, riskFreeRate, volatility);
            }
            sw.Stop();
            Console.WriteLine($"2000 binomial pricings:      {sw.Elapsed.TotalMilliseconds:F2} ms");
            Console.WriteLine($"Avg per pricing:             {sw.Elapsed.TotalMilliseconds / 2000:F4} ms\n");

            // Benchmark 4: Simulated Real-Time Feed
            Console.WriteLine("Benchmark 4: Simulated Real-Time Feed (1000 options/sec)");
            Console.WriteLine(new string('-', 60));

            int feedIterations = 10000;
            sw.Restart();
            for (int i = 0; i < feedIterations; i++)
            {
                // Simulate price update - calculate all Greeks
                var allGreeks = OptionPricing.CalculateAllGreeks(
                    underlyingPrice + (i % 10) * 0.01, // Slight price variation
                    optionType,
                    strikePrice,
                    timeToExpiration,
                    volatility,
                    riskFreeRate,
                    dividendYield);
            }
            sw.Stop();
            double feedTime = sw.Elapsed.TotalMilliseconds;
            double updatesPerSecond = feedIterations / (feedTime / 1000.0);
            Console.WriteLine($"Total time for {feedIterations} updates: {feedTime:F2} ms");
            Console.WriteLine($"Updates per second:          {updatesPerSecond:F0}");
            Console.WriteLine($"Avg time per update:         {feedTime / feedIterations:F4} ms\n");

            // Verification: Show actual values
            Console.WriteLine("Sample Calculation Results:");
            Console.WriteLine(new string('-', 60));
            var sampleGreeks = OptionPricing.CalculateAllGreeks(underlyingPrice, optionType, strikePrice,
                timeToExpiration, volatility, riskFreeRate, dividendYield);

            Console.WriteLine($"Option Price: ${sampleGreeks.Price:F4}");
            Console.WriteLine($"\nFirst Order Greeks:");
            Console.WriteLine($"  Delta:  {sampleGreeks.FirstOrder.Delta:F6}");
            Console.WriteLine($"  Gamma:  {sampleGreeks.FirstOrder.Gamma:F6}");
            Console.WriteLine($"  Theta:  {sampleGreeks.FirstOrder.Theta:F6}");
            Console.WriteLine($"  Vega:   {sampleGreeks.FirstOrder.Vega:F6}");
            Console.WriteLine($"  Rho:    {sampleGreeks.FirstOrder.Rho:F6}");
            Console.WriteLine($"\nSecond Order Greeks:");
            Console.WriteLine($"  Vanna:  {sampleGreeks.SecondOrder.Vanna:F6}");
            Console.WriteLine($"  Charm:  {sampleGreeks.SecondOrder.Charm:F6}");
            Console.WriteLine($"  Vomma:  {sampleGreeks.SecondOrder.Vomma:F6}");
            Console.WriteLine($"  Veta:   {sampleGreeks.SecondOrder.Veta:F6}");
            Console.WriteLine($"  Speed:  {sampleGreeks.SecondOrder.Speed:F6}");
            Console.WriteLine($"  Zomma:  {sampleGreeks.SecondOrder.Zomma:F6}");
            Console.WriteLine($"  Color:  {sampleGreeks.SecondOrder.Color:F6}");
        }
    }
}
