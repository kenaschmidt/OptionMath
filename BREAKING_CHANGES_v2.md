# Breaking Changes - Version 2.0

## Day Count Convention Accuracy Fix

### Summary

**BREAKING CHANGES** implemented to ensure accurate Greek calculations consistent with TradingCalendar's day count conventions. This fixes a **31% error in Theta** and other time-based Greeks for equity options.

### Date: February 15, 2026

---

## What Changed

### 1. Default Day Count Convention
```csharp
// BEFORE (v1.x):
public static int DayCountStandard { get; set; } = 365;

// AFTER (v2.0):
public static int DayCountStandard { get; set; } = 252;  // ⚠️ BREAKING CHANGE
```

**Impact:** Default changed from 365 to 252 to match equity options standard (most common use case).

### 2. Theta - Now Requires dayCountConvention Parameter

```csharp
// BEFORE (v1.x):
public static double Theta(
    double underlyingPrice, OptionType optionType, double optionStrikePrice,
    double timeToExpiration, double volatilityPercent,
    double riskFreeRatePercent, double dividendYieldPercent)

// AFTER (v2.0):
public static double Theta(
    double underlyingPrice, OptionType optionType, double optionStrikePrice,
    double timeToExpiration, double volatilityPercent,
    double riskFreeRatePercent, double dividendYieldPercent,
    int dayCountConvention = 252)  // ⚠️ NEW PARAMETER (with default)
```

**Impact:**
- Overloads with `DateTime + TradingInstrument` automatically use correct convention (no breaking change)
- Overloads with just `timeToExpiration` now have optional `dayCountConvention` parameter (defaults to 252)

### 3. Charm - Now Requires dayCountConvention Parameter

```csharp
// BEFORE (v1.x):
public static double Charm(
    double underlyingPrice, OptionType optionType, double optionStrikePrice,
    double timeToExpiration, double volatilityPercent,
    double riskFreeRatePercent, double dividendYieldPercent)

// AFTER (v2.0):
public static double Charm(
    double underlyingPrice, OptionType optionType, double optionStrikePrice,
    double timeToExpiration, double volatilityPercent,
    double riskFreeRatePercent, double dividendYieldPercent,
    int dayCountConvention = 252)  // ⚠️ NEW PARAMETER (with default)
```

### 4. Veta - Now Requires dayCountConvention Parameter

```csharp
// BEFORE (v1.x):
public static double Veta(
    double underlyingPrice, double optionStrikePrice,
    double timeToExpiration, double volatilityPercent,
    double riskFreeRatePercent, double dividendYieldPercent)

// AFTER (v2.0):
public static double Veta(
    double underlyingPrice, double optionStrikePrice,
    double timeToExpiration, double volatilityPercent,
    double riskFreeRatePercent, double dividendYieldPercent,
    int dayCountConvention = 252)  // ⚠️ NEW PARAMETER (with default)
```

### 5. Color - Now Requires dayCountConvention Parameter

```csharp
// BEFORE (v1.x):
public static double Color(
    double underlyingPrice, double optionStrikePrice,
    double timeToExpiration, double volatilityPercent,
    double riskFreeRatePercent, double dividendYieldPercent)

// AFTER (v2.0):
public static double Color(
    double underlyingPrice, double optionStrikePrice,
    double timeToExpiration, double volatilityPercent,
    double riskFreeRatePercent, double dividendYieldPercent,
    int dayCountConvention = 252)  // ⚠️ NEW PARAMETER (with default)
```

### 6. Batch Greeks Methods - Now Require dayCountConvention

```csharp
// BEFORE (v1.x):
public static FirstOrderGreeks CalculateFirstOrderGreeks(...)
public static SecondOrderGreeks CalculateSecondOrderGreeks(...)
public static CompleteGreeks CalculateAllGreeks(...)

// AFTER (v2.0):
public static FirstOrderGreeks CalculateFirstOrderGreeks(..., int dayCountConvention = 252)
public static SecondOrderGreeks CalculateSecondOrderGreeks(..., int dayCountConvention = 252)
public static CompleteGreeks CalculateAllGreeks(..., int dayCountConvention = 252)
```

---

## Why This Was Necessary

### The Problem

**Before v2.0:**
- TradingCalendar used 252 trading days/year for time-to-expiration calculation
- But Theta was normalized using 365 days/year
- Result: **31% error in Theta for equity options!**

**Example:**
```csharp
// Equity option with:
// - True theta per trading day: -$0.10
// - Old calculation (v1.x): -$0.069  (31% underestimation!)
// - New calculation (v2.0): -$0.10   (accurate!)

// Why? 365 / 252 = 1.448 error factor
```

### Impact by Asset Class

| Asset Class | Correct Days/Year | Old Error (v1.x) | Fixed in v2.0 |
|-------------|-------------------|------------------|---------------|
| **US Equities** | 252 | **-31% error** ❌ | ✅ Accurate |
| **EU Equities** | 250-256 | **-30 to -31% error** ❌ | ✅ Accurate |
| **FX Options** | 365 | ✅ Correct | ✅ Still correct |
| **Interest Rates** | 360 | **-1.4% error** ❌ | ✅ Accurate |

---

## Migration Guide

### If You Use TradingInstrument Overloads (✅ No Action Required)

```csharp
// This code continues to work WITHOUT changes:
var theta = OptionPricing.Theta(
    underlyingPrice, optionType, strikePrice,
    currentTime, expirationDate, instrument,  // ← TradingInstrument
    volatility, riskFreeRate, dividendYield);

// The correct day count is automatically inferred from the instrument!
```

### If You Use timeToExpiration Overloads (⚠️ Review Required)

#### Option 1: Accept New Default (Recommended for Equities)

```csharp
// No code changes needed - defaults to 252 (equity standard)
var theta = OptionPricing.Theta(
    underlyingPrice, optionType, strikePrice,
    timeToExpiration, volatility, riskFreeRate, dividendYield);
// ✅ Now uses 252 days/year (accurate for equities)
```

**Impact:** Your Theta values will **increase by ~31%** (now correct!)

#### Option 2: Explicitly Specify Convention

```csharp
// For FX options (calendar days):
var theta = OptionPricing.Theta(
    underlyingPrice, optionType, strikePrice,
    timeToExpiration, volatility, riskFreeRate, dividendYield,
    dayCountConvention: 365);  // ← Explicit

// For bonds/money market (360 days):
var theta = OptionPricing.Theta(
    underlyingPrice, optionType, strikePrice,
    timeToExpiration, volatility, riskFreeRate, dividendYield,
    dayCountConvention: 360);  // ← Explicit
```

#### Option 3: Maintain Old Behavior (Not Recommended)

```csharp
// To get old v1.x behavior (365 days):
var theta = OptionPricing.Theta(
    underlyingPrice, optionType, strikePrice,
    timeToExpiration, volatility, riskFreeRate, dividendYield,
    dayCountConvention: 365);  // ⚠️ Old behavior (less accurate for equities)
```

### Batch Calculations

```csharp
// BEFORE (v1.x):
var greeks = OptionPricing.CalculateAllGreeks(
    underlyingPrice, optionType, strikePrice,
    timeToExpiration, volatility, riskFreeRate, dividendYield);

// AFTER (v2.0) - Same call works, but uses 252 default:
var greeks = OptionPricing.CalculateAllGreeks(
    underlyingPrice, optionType, strikePrice,
    timeToExpiration, volatility, riskFreeRate, dividendYield);
// ✅ Now uses 252 days/year

// Or explicit for FX:
var greeks = OptionPricing.CalculateAllGreeks(
    underlyingPrice, optionType, strikePrice,
    timeToExpiration, volatility, riskFreeRate, dividendYield,
    dayCountConvention: 365);  // ← For FX options
```

---

## What's NOT Affected

These Greeks **don't use day count** and are unchanged:

✅ **Delta** - No change
✅ **Gamma** - No change
✅ **Vega** - No change
✅ **Rho** - No change
✅ **Vanna** - No change
✅ **Vomma** - No change
✅ **Zomma** - No change
✅ **Speed** - No change
✅ **Option Pricing** - No change

Only time-based Greeks affected:
- ⚠️ Theta
- ⚠️ Charm
- ⚠️ Veta
- ⚠️ Color

---

## New Helper Method: GetDayCountConvention

The library now automatically determines day count from instrument type:

```csharp
// Internal method (you don't call this, but it's what happens):
private static int GetDayCountConvention(TradingInstrument instrument)
{
    // FX options: 365 calendar days
    if (instrument.Name.Contains("fx")) return 365;

    // Bonds/rates: 360 days (money market)
    if (instrument.Name.Contains("rate") || instrument.Name.Contains("bond")) return 360;

    // Default: Equities use 252 trading days
    return 252;
}
```

---

## Expected Changes in Your Results

### For Equity Options (Most Common):

| Greek | Old Value (v1.x) | New Value (v2.0) | Change |
|-------|------------------|------------------|--------|
| **Theta** | -0.069 | -0.100 | **+45% (more negative)** ✅ |
| **Charm** | -0.0002 | -0.0003 | **+45%** ✅ |
| **Veta** | -0.0008 | -0.0012 | **+45%** ✅ |
| **Color** | 0.00012 | 0.00017 | **+45%** ✅ |

**Why +45%?** Because 365/252 = 1.448, so the denominator is smaller, making the absolute values larger.

### For FX Options (Using dayCountConvention: 365):

✅ **No change** - Results identical to v1.x

---

## Testing Your Migration

### Step 1: Update Your Code

Add `dayCountConvention` parameters where needed (or rely on defaults).

### Step 2: Validate Results

```csharp
// Test case - Equity option
double theta_v2 = OptionPricing.Theta(..., dayCountConvention: 252);
double theta_old = OptionPricing.Theta(..., dayCountConvention: 365);

// Expected relationship:
// theta_v2 ≈ theta_old * (365 / 252) = theta_old * 1.448

Console.WriteLine($"Old: {theta_old}, New: {theta_v2}, Ratio: {theta_v2/theta_old}");
// Should print ratio ≈ 1.45
```

### Step 3: Update Tests

If you have hardcoded expected values in tests:

```csharp
// BEFORE (v1.x test):
Assert.AreEqual(-0.069, theta, 0.001);

// AFTER (v2.0 test) - for equities:
Assert.AreEqual(-0.100, theta, 0.001);  // ← Updated expected value
```

---

## Frequently Asked Questions

### Q: Will my existing code break?

**A:** Most code will compile without changes due to default parameters. However:
- ✅ Code using `DateTime + TradingInstrument` overloads: **No changes needed**
- ⚠️ Code using `timeToExpiration` overloads: **Will compile, but results change**
- ⚠️ Hardcoded test assertions: **May need updating**

### Q: Should I use 252 or 365?

**A:**
- **252** = Equities, index options, ETF options (US/most markets)
- **365** = FX options, currency pairs
- **360** = Money market instruments, some bonds
- **250-256** = Equities in some European markets

### Q: How do I know what my current code was using?

**A:** In v1.x, everything used 365. So if you were trading equities, your Theta was **31% too small**.

### Q: What if I have mixed asset types?

**A:** Use the `TradingInstrument` overloads - they automatically use the correct convention for each instrument type.

### Q: Can I revert to old behavior?

**A:** Yes, but not recommended:
```csharp
OptionPricing.DayCountStandard = 365;  // Global change (affects all calls)
// OR
Theta(..., dayCountConvention: 365);   // Per-call basis
```

---

## Benefits of v2.0

✅ **Accurate Greeks** for equity options (31% error eliminated)
✅ **Consistent with TradingCalendar** approach
✅ **Instrument-aware** - automatically uses correct convention
✅ **Flexible** - can override for special cases
✅ **Production-ready** for professional trading systems

---

## Rollout Recommendation

### Phase 1: Test Environment
1. Update to v2.0 in test/dev environment
2. Compare new vs old Theta values
3. Validate that new values are ~45% larger for equities
4. Update any hardcoded test expectations

### Phase 2: Validation
1. Run side-by-side with production for 1 week
2. Monitor that equity Theta is now properly valued
3. Verify FX options (if using 365) produce same results

### Phase 3: Production
1. Deploy v2.0 to production
2. Monitor P&L attribution
3. Confirm time decay matches expected behavior

---

## Support

If you encounter issues during migration:

1. Check `DAY_COUNT_ANALYSIS.md` for detailed explanation
2. Review your asset types and select appropriate convention
3. Use the migration examples above
4. Run comparative tests between v1.x and v2.0

---

**Version:** 2.0.0
**Breaking Change Severity:** HIGH (affects accuracy of time-based Greeks)
**Recommended Action:** Update and test thoroughly before production deployment
