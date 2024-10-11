# {py:mod}`gauss_markov_srs.gauss_markov`

```{py:module} gauss_markov_srs.gauss_markov
```

```{autodoc2-docstring} gauss_markov_srs.gauss_markov
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`FitResult <gauss_markov_srs.gauss_markov.FitResult>`
  - ```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitResult
    :summary:
    ```
* - {py:obj}`FitStats <gauss_markov_srs.gauss_markov.FitStats>`
  - ```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`gauss_markov_fit <gauss_markov_srs.gauss_markov.gauss_markov_fit>`
  - ```{autodoc2-docstring} gauss_markov_srs.gauss_markov.gauss_markov_fit
    :summary:
    ```
* - {py:obj}`fit_stats <gauss_markov_srs.gauss_markov.fit_stats>`
  - ```{autodoc2-docstring} gauss_markov_srs.gauss_markov.fit_stats
    :summary:
    ```
````

### API

`````{py:class} FitResult
:canonical: gauss_markov_srs.gauss_markov.FitResult

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitResult
```

````{py:attribute} q
:canonical: gauss_markov_srs.gauss_markov.FitResult.q
:type: numpy.ndarray
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitResult.q
```

````

````{py:attribute} Cqq
:canonical: gauss_markov_srs.gauss_markov.FitResult.Cqq
:type: numpy.ndarray
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitResult.Cqq
```

````

````{py:attribute} qu
:canonical: gauss_markov_srs.gauss_markov.FitResult.qu
:type: numpy.ndarray
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitResult.qu
```

````

````{py:attribute} weights
:canonical: gauss_markov_srs.gauss_markov.FitResult.weights
:type: numpy.ndarray
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitResult.weights
```

````

`````

`````{py:class} FitStats
:canonical: gauss_markov_srs.gauss_markov.FitStats

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats
```

````{py:attribute} n_meas
:canonical: gauss_markov_srs.gauss_markov.FitStats.n_meas
:type: int
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.n_meas
```

````

````{py:attribute} n_adj
:canonical: gauss_markov_srs.gauss_markov.FitStats.n_adj
:type: int
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.n_adj
```

````

````{py:attribute} n_excluded
:canonical: gauss_markov_srs.gauss_markov.FitStats.n_excluded
:type: int
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.n_excluded
```

````

````{py:attribute} n_included
:canonical: gauss_markov_srs.gauss_markov.FitStats.n_included
:type: int
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.n_included
```

````

````{py:attribute} n_dof
:canonical: gauss_markov_srs.gauss_markov.FitStats.n_dof
:type: int
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.n_dof
```

````

````{py:attribute} y
:canonical: gauss_markov_srs.gauss_markov.FitStats.y
:type: numpy.ndarray
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.y
```

````

````{py:attribute} yu
:canonical: gauss_markov_srs.gauss_markov.FitStats.yu
:type: numpy.ndarray
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.yu
```

````

````{py:attribute} self_sensitivity
:canonical: gauss_markov_srs.gauss_markov.FitStats.self_sensitivity
:type: numpy.ndarray
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.self_sensitivity
```

````

````{py:attribute} residuals
:canonical: gauss_markov_srs.gauss_markov.FitStats.residuals
:type: numpy.ndarray
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.residuals
```

````

````{py:attribute} relative_residuals
:canonical: gauss_markov_srs.gauss_markov.FitStats.relative_residuals
:type: numpy.ndarray
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.relative_residuals
```

````

````{py:attribute} chi2
:canonical: gauss_markov_srs.gauss_markov.FitStats.chi2
:type: float
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.chi2
```

````

````{py:attribute} birge
:canonical: gauss_markov_srs.gauss_markov.FitStats.birge
:type: float
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.birge
```

````

````{py:attribute} birge_limit
:canonical: gauss_markov_srs.gauss_markov.FitStats.birge_limit
:type: float
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.birge_limit
```

````

````{py:attribute} p_value
:canonical: gauss_markov_srs.gauss_markov.FitStats.p_value
:type: float
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.p_value
```

````

````{py:attribute} n_corr_nonzero
:canonical: gauss_markov_srs.gauss_markov.FitStats.n_corr_nonzero
:type: int
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.n_corr_nonzero
```

````

````{py:attribute} min_eigenvalue
:canonical: gauss_markov_srs.gauss_markov.FitStats.min_eigenvalue
:type: float
:value: >
   None

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.FitStats.min_eigenvalue
```

````

`````

````{py:function} gauss_markov_fit(A, y, Cyy, mask=None)
:canonical: gauss_markov_srs.gauss_markov.gauss_markov_fit

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.gauss_markov_fit
```
````

````{py:function} fit_stats(A, y, Cyy, fit_result, mask=None)
:canonical: gauss_markov_srs.gauss_markov.fit_stats

```{autodoc2-docstring} gauss_markov_srs.gauss_markov.fit_stats
```
````
