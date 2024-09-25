import gauss_markov_srs as gm

input_data = gm.load_excel("./Data/Input/CCTF2021/inputs.xlsx", dcols=[5, 6])
input_corr = gm.load_txt("./Data/Input/CCTF2021/correlations.txt")
reference = gm.load_txt("./Data/Reference/CCTF2017.txt", dcols=[2])


input_data = gm.data_with_long_id(input_data)


y, yu = gm.get_y(input_data, reference)
M = gm.get_model_matrix(input_data, reference)
cyy = gm.get_corr_matrix(input_data, input_corr)

Cyy = gm.corr2cov(cyy, yu)


fit_results = gm.gauss_markov_fit(M, y, Cyy)
stats = gm.fit_stats(M, y, Cyy, fit_results)


gm.pretty_print_fit_result(fit_results, reference)
gm.pretty_print_stats(stats)

Sc = stats.self_sensitivity

mask = Sc > 0.01


fit_results2 = gm.gauss_markov_fit(M, y, Cyy, mask=mask)
stats2 = gm.fit_stats(M, y, Cyy, fit_results2, mask=mask)

print("\n====\n")
gm.pretty_print_fit_result(fit_results2, reference)
gm.pretty_print_stats(stats2)
