# Employment distribution in EUROSTAT countries

According to the three-sector theory, employment shifts from the primary
sector (raw material extraction), to the secondary sector (industry,
energy, and construction), and then to the tertiary sector (services) as
economies develop. The \`eurostat_employment\` data set contains
EUROSTAT data on employment, aggregated for both sexes and all ages,
distributed by economic activity in 2008 for 29 EUROSTAT member
countries.

A related variable is the logarithm of gross domestic product per person
in EUR at current prices (\`logGDP\`). For exploratory purposes, it is
also categorised as a binary variable indicating values above or below
the median (\`Binary GDP\`).

The employment composition has 11 parts:

- Primary sector

- Manufacturing

- Energy

- Construction

- Trade repair transport

- Hotels restaurants

- Financial intermediation

- Real estate

- Educ admin defense soc sec

- Health social work

- Other services

## Usage

``` r
eurostat_employment
```

## Format

An object of class `data.frame` with 29 rows and 17 columns.
