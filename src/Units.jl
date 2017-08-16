using Unitful
using Unitful.DefaultSymbols
Unitful.register(current_module());

@unit day "day" Day (24 * u"hr") false
@unit week "week" Week (7 * u"day") false
@unit month "month" Month ((365/12) * u"day") false
@unit year "year" Year (12 * u"month") false
@unit year "year" Year (365 * u"day") false

Rates = Union{
  Unitful.ContextUnits{N,Unitful.Dimensions{(Unitful.Dimension{:Time}(-1//1),)},P} where P where N,
  Unitful.FixedUnits{N,Unitful.Dimensions{(Unitful.Dimension{:Time}(-1//1),)}} where N,
  Unitful.FreeUnits{N,Unitful.Dimensions{(Unitful.Dimension{:Time}(-1//1),)}} where N}
