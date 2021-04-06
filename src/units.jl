module Units

import Unitful
using Unitful: @unit
day = Unitful.d
week = Unitful.wk
year = Unitful.yr
@unit month "month" Month 2.628e6 * Unitful.s false

const days = day
const weeks = week
const months = month
const years = year
const January = 0month
const February = 1month
const March = 2months
const April = 3months
const May = 4months
const June = 5months
const July = 6months
const August = 7months
const September = 8months
const October = 9months
const November = 10months
const December = 11months

const localunits = Unitful.basefactors
function __init__()
    merge!(Unitful.basefactors, localunits)
    Unitful.register(Units)
end

export day, days, week, weeks, month, months, year, years, Rates,
January, February, March, April, May, June, July, August,
September, October, November, December

end
