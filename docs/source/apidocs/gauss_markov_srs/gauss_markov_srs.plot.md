# {py:mod}`gauss_markov_srs.plot`

```{py:module} gauss_markov_srs.plot
```

```{autodoc2-docstring} gauss_markov_srs.plot
:allowtitles:
```

## Module Contents

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`plot_correlation_matrix <gauss_markov_srs.plot.plot_correlation_matrix>`
  - ```{autodoc2-docstring} gauss_markov_srs.plot.plot_correlation_matrix
    :summary:
    ```
* - {py:obj}`plot_residuals <gauss_markov_srs.plot.plot_residuals>`
  - ```{autodoc2-docstring} gauss_markov_srs.plot.plot_residuals
    :summary:
    ```
* - {py:obj}`plot_fit_results <gauss_markov_srs.plot.plot_fit_results>`
  - ```{autodoc2-docstring} gauss_markov_srs.plot.plot_fit_results
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`params <gauss_markov_srs.plot.params>`
  - ```{autodoc2-docstring} gauss_markov_srs.plot.params
    :summary:
    ```
````

### API

````{py:data} params
:canonical: gauss_markov_srs.plot.params
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.plot.params
```

````

````{py:function} plot_correlation_matrix(dir, label, cyy, figname='correlation.png')
:canonical: gauss_markov_srs.plot.plot_correlation_matrix

```{autodoc2-docstring} gauss_markov_srs.plot.plot_correlation_matrix
```
````

````{py:function} plot_residuals(dir, label, data, reference, fit_result: gauss_markov_srs.gauss_markov.FitResult, fit_stats: gauss_markov_srs.gauss_markov.FitStats, mask=None, basename='residuals{}-{}.png')
:canonical: gauss_markov_srs.plot.plot_residuals

```{autodoc2-docstring} gauss_markov_srs.plot.plot_residuals
```
````

````{py:function} plot_fit_results(dir, reference, fit_results, labels, figname='fit-result.png')
:canonical: gauss_markov_srs.plot.plot_fit_results

```{autodoc2-docstring} gauss_markov_srs.plot.plot_fit_results
```
````
