using .PyCall
py"""
from math import floor
from ecmwfapi import ECMWFDataServer

server = ECMWFDataServer()


def retrieve_ECMWF(param, from_year, to_year, filename, **kwargs):

    months = range(1, 13)
    years = range(from_year, to_year + 1)
    decades = sorted({floor(y / 10) * 10 for y in years})

    # Loop through decades and create a request list for all months/years
    for d in decades:
        # Filter for years within the decade
        years_in_decade = list(filter(lambda y: floor(y / 10) * 10 == d, years))

        # Set up all request months per decade in correct format
        request_dates = "/".join([f'{y}{m:02}01' for y in years_in_decade for m in months])

        # Create target file
        target = f'{filename}_{d}'
        print(f'Years: {years_in_decade}\nOutput file: {target}')
        ECMWF_request(param, request_dates, d, target, **kwargs)


def ECMWF_request(
        param, request_dates, decade, target,
        eclass = 'ei', dataset = 'interim', stream='moda', modeltype='an', levtype='sfc', grid='0.75/0.75', format='netcdf', step='0-12', num = '0'
):

    server.retrieve({
        'class':   eclass,
        'dataset': dataset,
        'expver':  '1',
        'stream':  stream,
        'number':  num,
        'type':    modeltype,
        'levtype': levtype,
        'param':   param,
        'grid':    grid,
        'format':  format,
        'date':    request_dates,
        'decade':  decade,
        'target':  target,
        'step':    step
    })
"""

function retrieve_ECMWF(param::String, from_year::Int64, to_year::Int64, filename::String = "era_interim"; kws...)
    py"retrieve_ECMWF"(param, from_year, to_year, filename; kws...)
end
