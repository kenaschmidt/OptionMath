# Day Count Convention Inconsistency Analysis

## Issue Identified

The OptionMath library has an inconsistency between:
1. **TradingCalendar's sophisticated time calculations** (used for time-to-expiration)
2. **DayCountStandard = 365** (used for normalizing Greeks to "per day")

## Where DayCountStandard is Used

### 1. Theta Normalization (Lines 616, 628, 648, 657)
```csharp
// Theta formula gives change per year
// We divide by DayCountStandard to get "per day"
return ret / DayCountStandard;  // Uses 365
```

### 2. Time-Based Greeks: Charm, Veta, Color (Lines 870, 879, 966, 973, 1108, 1117)
```csharp
// Check if less than 1 day remaining
if (timeToExpiration <= (1.0 / DayCountStandard))

// Calculate 1-day difference
double delta2 = Delta(..., timeToExpiration - (1.0 / DayCountStandard), ...);
```

## The Problem

### Example Scenario:
**Equity Options** (common convention: 252 trading days/year)

1. **Time to Expiration Calculation:**
   - TradingCalendar uses **252 trading days** per year
   - 30 calendar days with 20 trading days = 20/252 = 0.0794 years

2. **Theta "Per Day" Calculation:**
   - Theta formula calculates change per year
   - We divide by **365** to get "per day"
   - But if the year has 252 trading days, we should divide by 252!

### Impact on Theta:
```
Correct Theta (per trading day) = Theta_annual / 252
Current Theta (per calendar day) = Theta_annual / 365

Error Factor = 365 / 252 = 1.448

Example:
- True theta per trading day: -$0.10
- Reported theta (current): -$0.069
- **31% underestimation!**
```

### Impact on Charm, Veta, Color:
These use numerical differentiation with a 1-day step:
```csharp
timeToExpiration - (1.0 / DayCountStandard)
// Subtracts 1/365 year instead of 1/252 year
// Results in ~31% error in the time step
```

## Industry Conventions

Different markets use different conventions:

| Asset Class | Day Count Convention | Days Per Year |
|-------------|---------------------|---------------|
| **Equity Options** | Trading days (US) | 252 |
| **Equity Options** | Trading days (Europe) | 250-256 |
| **Index Options** | Trading days | 252 |
| **FX Options** | Calendar days | 365 |
| **Interest Rate Options** | Actual/360 or Actual/365 | 360 or 365 |
| **Commodity Options** | Varies by commodity | 252-365 |

## Current Behavior Analysis

### What Works Correctly:
✅ **Time to expiration** - Uses TradingCalendar's sophisticated calculation
✅ **Option pricing** - Correctly uses time to expiration
✅ **Delta, Gamma, Vega, Rho** - Not affected (don't use DayCountStandard)

### What's Potentially Incorrect:
⚠️ **Theta** - Normalized using 365, should match TradingCalendar's convention
⚠️ **Charm** - Uses 1/365 year step, should use convention-appropriate step
⚠️ **Veta** - Uses 1/365 year step, should use convention-appropriate step
⚠️ **Color** - Uses 1/365 year step, should use convention-appropriate step

## Severity Assessment

### High Impact Scenarios:
1. **Trading Systems** that use Theta for hedging decisions
2. **Risk Management** systems that aggregate Theta across portfolios
3. **P&L Attribution** that decomposes time decay
4. **Automated Market Making** that relies on accurate Greeks

### Low Impact Scenarios:
1. **Directional Trading** (Delta is correct)
2. **Volatility Trading** (Vega is correct)
3. **Spread Trading** (if both legs use same convention, error cancels)

## Solutions

### Option 1: Make DayCountStandard Configurable Per Instrument ✓ RECOMMENDED
```csharp
public static double Theta(..., int? dayCountConvention = null)
{
    int dayCount = dayCountConvention ?? DayCountStandard;
    return ret / dayCount;
}
```

**Pros:**
- Backward compatible (default to 365)
- Allows precision when needed
- User controls the convention

**Cons:**
- Requires user knowledge of correct convention
- More parameters to manage

### Option 2: Integrate with TradingCalendar's Convention
```csharp
// Add to TradingInstrument class:
public int DayCountConvention { get; set; }

// Use in calculations:
public static double Theta(..., TradingCalendar.TradingInstrument instrument)
{
    int dayCount = instrument.DayCountConvention;
    return ret / dayCount;
}
```

**Pros:**
- Fully consistent with TradingCalendar
- Instrument-specific conventions
- Single source of truth

**Cons:**
- Requires TradingCalendar changes
- All overloads need instrument parameter
- More breaking changes

### Option 3: Add Explicit "PerTradingDay" vs "PerCalendarDay" Methods
```csharp
public static double ThetaPerCalendarDay(...)
{
    return CalculateTheta(...) / 365;
}

public static double ThetaPerTradingDay(..., int tradingDaysPerYear = 252)
{
    return CalculateTheta(...) / tradingDaysPerYear;
}
```

**Pros:**
- Explicit about what you're getting
- No ambiguity
- User chooses

**Cons:**
- Doubles the number of methods
- Existing code unclear about which it uses

### Option 4: Document Current Behavior, Add Warning ⚠️ MINIMUM
```csharp
/// <summary>
/// Calculates Theta - change in option price per calendar day.
/// NOTE: Result is normalized using 365 calendar days per year.
/// For trading days, multiply by (365 / tradingDaysPerYear).
/// </summary>
public static double Theta(...)
```

**Pros:**
- No code changes
- Backward compatible
- Users can adjust if needed

**Cons:**
- Doesn't fix the problem
- Burden on users

## Recommended Approach

### Phase 1: Immediate (Backward Compatible)
1. **Document the current behavior** clearly in XML comments
2. **Add optional dayCount parameter** to Theta, Charm, Veta, Color:
   ```csharp
   public static double Theta(..., int dayCountConvention = 365)
   ```
3. **Add helper methods** for common conventions:
   ```csharp
   public static double ThetaTradingDays(...) => Theta(..., 252);
   ```

### Phase 2: Future Enhancement
1. **Extend TradingInstrument** to include day count convention
2. **Add overloads** that automatically use instrument's convention
3. **Deprecate** (but don't remove) the fixed 365 approach

## Testing Requirements

If we fix this, we need to:
1. ✅ Verify Theta matches expected values for 252-day convention
2. ✅ Test that 365-day convention still works for FX/rates
3. ✅ Ensure Charm, Veta, Color use consistent day steps
4. ✅ Add tests comparing different conventions
5. ✅ Document which convention each test uses

## Example Fix (Option 1)

```csharp
/// <summary>
/// Calculates Theta - change in option price per day
/// </summary>
/// <param name="dayCountConvention">Days per year (252 for equities, 365 for FX, etc.)</param>
public static double Theta(
    double underlyingPrice,
    OptionType optionType,
    double optionStrikePrice,
    double timeToExpiration,
    double volatilityPercent,
    double riskFreeRatePercent,
    double dividendYieldPercent,
    int dayCountConvention = 365)  // Default maintains backward compatibility
{
    // ... existing calculation ...

    return ret / dayCountConvention;  // Use provided convention
}

// Convenience methods
public static double ThetaEquities(...) => Theta(..., dayCountConvention: 252);
public static double ThetaFX(...) => Theta(..., dayCountConvention: 365);
```

## Impact on Existing Code

### If We Add Optional Parameter (Recommended):
- ✅ **Existing code continues to work** (uses default 365)
- ✅ **New code can specify convention** for accuracy
- ✅ **No breaking changes**
- ⚠️ **Existing code may be using wrong convention** (but at least consistently)

### If We Change Default to 252:
- ❌ **Breaking change** - all existing Theta values change by 31%
- ❌ **Existing tests would fail**
- ❌ **User code would break**

## Conclusion

**YES, this is affecting calculations and could create problems.**

### Severity:
- **Medium to High** for trading/risk systems using Theta
- **Low** for systems only using Delta, Gamma, Vega (unaffected)

### Recommendation:
1. **Short term:** Add optional `dayCountConvention` parameter (default 365 for compatibility)
2. **Medium term:** Document which convention to use for each asset class
3. **Long term:** Integrate with TradingCalendar's instrument-specific conventions

### Action Items:
1. ✅ Document current behavior clearly
2. ✅ Add optional dayCountConvention parameter to affected methods
3. ✅ Add convenience methods for common conventions (252, 365, 360)
4. ✅ Update tests to validate both conventions
5. ✅ Add warnings to PERFORMANCE_OPTIMIZATIONS.md

---

**Status:** Issue Confirmed - Solution Proposed
**Priority:** Medium-High (affects accuracy for equity options traders)
**Backward Compatibility:** Maintained with optional parameters
