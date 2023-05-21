import pandas as pd
import numpy as np
import utils
from tqdm import tqdm
from time import time


def gen_nans_mask(n, m, p):
    """
    @params m, n: width and length of the mask
    @param p: nans fraction
    @return: nans mask
    """
    nans_mask = np.random.binomial(1, p, size=(n, m - 2))
    nans_mask = np.concatenate((np.zeros((n, 2)), nans_mask), axis=1).astype(bool)
    return nans_mask


def add_metrics(diffs, rmse_lst, mae_lst):
    """
    @param diffs: differences
    @param rmse_lst: list accumulating all rmse
    @param mae_lst: list accumulating all mae
    """
    rmse_lst.append(np.sqrt((diffs ** 2).mean()))
    mae_lst.append(np.abs(diffs).mean())


def methyLImp_eps_dependence(filename, eps_arr, p=0.01, p_impute=1e-2, verbose=1):
    """
    @param filename: path to the file with data
    @param eps_arr: array of half neighborhood sisez that need to be tested
    @param p: nans fraction to test
    @param p_impute: fraction of nans that will be imputed
    @param verbose: verbosity of the function
    @return: lists with rmse, mae, and time
    """
    data  = pd.read_csv(filename, sep='\t')
    mask = np.invert(np.any(np.isnan(data.values[:, 2:].astype(np.float64)), axis=1))
    data = data[mask].reset_index(drop=True)
    values = data.values

    rmse_lst = []
    mae_lst = []
    time_lst = []

    if verbose:
        eps_arr = tqdm(eps_arr)
    for eps in eps_arr:
        nans_mask = gen_nans_mask(*data.shape, p)
        data_with_nans = data.copy()
        data_with_nans[nans_mask] = np.nan

        rows = np.random.binomial(1, p_impute, size=data.shape[0]).astype(bool)
        rows = rows & np.any(nans_mask, axis=1)
        impute_nans_mask = nans_mask & rows[:, None]
        
        t_start = time()
        imputed = utils.methyLImp(data_with_nans, eps=eps, rows=rows, verbose=max(verbose - 1, 0))
        time_lst.append(time() - t_start)

        diffs_methyLImp = imputed[impute_nans_mask[:, 2:]].astype(float) - values[impute_nans_mask].astype(float)
        add_metrics(diffs_methyLImp, rmse_lst, mae_lst)

    return rmse_lst, mae_lst, time_lst


def nbp_eps_dependence(filename, eps_arr, p=0.01, p_impute=1e-3, verbose=1):
    """
    @param filename: path to the file with data
    @param eps_arr: array of half neighborhood sisez that need to be tested
    @param p: nans fraction to test
    @param p_impute: fraction of nans that will be imputed
    @param verbose: verbosity of the function
    @return: lists with rmse, mae, and time
    """
    data  = pd.read_csv(filename, sep='\t')
    mask = np.invert(np.any(np.isnan(data.values[:, 2:].astype(np.float64)), axis=1))
    data = data[mask].reset_index(drop=True)
    values = data.values

    rmse_lst = []
    mae_lst = []
    time_lst = []

    if verbose:
        eps_arr = tqdm(eps_arr)
    for eps in eps_arr:
        nans_mask = gen_nans_mask(*data.shape, p)
        data_with_nans = data.copy()
        data_with_nans[nans_mask] = np.nan

        rows = np.random.binomial(1, p_impute, size=data.shape[0]).astype(bool)
        rows = rows & np.any(nans_mask, axis=1)
        impute_nans_mask = nans_mask & rows[:, None]
        
        t_start = time()
        imputed = utils.impute_nbp(data_with_nans, eps=eps, impute_positions=impute_nans_mask[:, 2:], verbose=max(verbose - 1, 0))
        time_lst.append(time() - t_start)

        diffs_methyLImp = imputed.values[impute_nans_mask].astype(float) - values[impute_nans_mask].astype(float)
        add_metrics(diffs_methyLImp, rmse_lst, mae_lst)

    return rmse_lst, mae_lst, time_lst


def create_benchmark_table(*args, **kwargs):
    """
    @params the same as for run_experiments
    """
    results = run_experiments(*args, **kwargs)
    for metric in list(results.keys())[1:]:
        results[metric] = [f"{round(np.mean(lst), 2)} +/- {round(np.std(lst), 2)}" for lst in results[metric]]
    benchmark_table = pd.DataFrame(results)
    return benchmark_table.drop(benchmark_table.columns[-4:], axis=1), benchmark_table.drop(benchmark_table.columns[1:5], axis=1)


def run_experiments(filename, N_iter, eps_methyLImp, eps_nbp, p_arr=[0.0005, 0.001, 0.005, 0.01], p_impute=1e-3, verbose=1, gen_nm=gen_nans_mask):
    """
    @param filename: path to the file with data
    @param N_iter: number of iterations for each set of hyperparameters
    @param eps_methyLImp: half neighborhood size for the methyLImp
    @param eps_nbp: half neighborhood size for the nbp
    @param p_arr: values of nans fraction
    @param p_impute: fraction of nans that will be imputed
    @param verbose: verbosity of the function
    @param gen_nm: function which generates nans mask
    @return: benchmark table with results of the experiments
    """
    data  = pd.read_csv(filename, sep='\t')
    mask = np.invert(np.any(np.isnan(data.values[:, 2:].astype(np.float64)), axis=1))
    data = data[mask].reset_index(drop=True)
    values = data.values

    benchmark_table = {
        "nans percent": [],
        "methyLImp rmse": [],
        "nbp rmse": [],
        "cytosine mean rmse": [],
        "people mean rmse": [],
        "methyLImp mae": [],
        "nbp mae": [],
        "cytosine mean mae": [],
        "people mean mae": [],
    }
    keys = benchmark_table.keys()
    
    if verbose:
        p_arr = tqdm(p_arr)
    for p in p_arr:
        metrics_lists = [[] for _ in range(8)]
        if verbose >= 2:
            lst = tqdm(range(N_iter))
        else:
            lst = range(N_iter)
        for _ in lst:
            nans_mask = gen_nm(*data.shape, p)
            data_with_nans = data.copy()
            data_with_nans[nans_mask] = np.nan

            rows = np.random.binomial(1, p_impute, size=data.shape[0]).astype(bool)
            rows = rows & np.any(nans_mask, axis=1)
            impute_nans_mask = nans_mask & rows[:, None]

            
            imputed_methyLImp = utils.methyLImp(data_with_nans, eps=eps_methyLImp, rows=rows, verbose=max(verbose - 2, 0))
            diffs_methyLImp = imputed_methyLImp[impute_nans_mask[:, 2:]].astype(float) - values[impute_nans_mask].astype(float)
            add_metrics(diffs_methyLImp, metrics_lists[0], metrics_lists[4])

            imputed_nbp = utils.impute_nbp(data_with_nans, eps=eps_nbp, impute_positions=impute_nans_mask[:, 2:], verbose=max(verbose - 2, 0))
            diffs_nbp = imputed_nbp.values[impute_nans_mask].astype(float) - values[impute_nans_mask].astype(float)
            add_metrics(diffs_nbp, metrics_lists[1], metrics_lists[5])

            diffs_cytosine_mean = (np.nanmean(data_with_nans.values[:, 2:].astype(float), axis=1)[:, None] - values[:, 2:].astype(float))[impute_nans_mask[:, 2:]]
            add_metrics(diffs_cytosine_mean, metrics_lists[2], metrics_lists[6])

            diffs_people_mean = (np.nanmean(data_with_nans.values[:, 2:].astype(float), axis=0)[None, :] - values[:, 2:].astype(float))[impute_nans_mask[:, 2:]]
            add_metrics(diffs_people_mean, metrics_lists[3], metrics_lists[7])
        
        for key, val in zip(keys, [p] + metrics_lists):
            benchmark_table[key].append(val)

    return benchmark_table


def combo_errors(filename, eps_methyLImp, eps_nbp1, eps_nbp2, N_iter=10, p=0.05):
    """
    @param filename: path to the file with data
    @param eps_methyLImp: half neighborhood size for the methyLImp
    @param eps_nbp1: half neighborhood size for the nbp in the first part of the algorithm
    @param eps_nbp2: half neighborhood size for the nbp in the second part of the algorithm
    @param N_iter: number of iterations for each set of hyperparameters
    @param p: nans fraction
    @return: lists with rmse and mae for all sets of hyperparameters
    """
    data  = pd.read_csv(filename, sep='\t')
    mask = np.invert(np.any(np.isnan(data.values[:, 2:].astype(np.float64)), axis=1))
    data_without_nans = data[mask].reset_index(drop=True)

    n, m = data_without_nans.shape
    
    rmse_lst = []
    mae_lst = []

    for _ in range(N_iter):
        row_inds = np.random.choice(np.arange(n), size=n // 5000, replace=False)
        data_without_rows = data_without_nans.drop(row_inds).reset_index(drop=True)

        n, m = data_without_rows.shape

        nans_mask = gen_nans_mask(n, m, p)
        data_with_nans = data_without_rows.copy()
        data_with_nans[nans_mask] = np.nan

        slices = []

        chromosome_coords = utils.get_chromosome_coords(data_with_nans.chromosome)
        positions = data_with_nans.position.values

        for row in row_inds:
            c = data_without_nans.chromosome[row]
            pos = data_without_nans.position[row]
            i, j = chromosome_coords[c]
            l, r = utils.find_slice(pos, eps_nbp2, positions[i:j])
            slices.append((l + i, r + i))

        rows_to_impute = np.zeros(len(data_with_nans), dtype=bool)
        for l, r in slices:
            rows_to_impute[l:r] = True

        # which missing values have to be imputed
        impute_nans_mask = nans_mask & rows_to_impute[:, None]
        print(impute_nans_mask.sum())

        print("methyLImp")
        imputed_data_methyLimp = utils.methyLImp(data_with_nans, eps=eps_methyLImp, rows=rows_to_impute & np.any(nans_mask, axis=1))
        
        print("nbp")
        imputed_data_nbp = utils.impute_nbp(data_with_nans, eps=eps_nbp1, impute_positions=impute_nans_mask[:, 2:])
        
        print("people mean")
        imputed_data_people_mean = utils.impute_people_mean(data_with_nans, impute_positions=impute_nans_mask[:, 2:])

        print("cytosine mean")
        imputed_data_cytosine_mean = utils.impute_cytosine_mean(data_with_nans, impute_positions=impute_nans_mask[:, 2:])
        

        def impute_rows(imputed_data):
            rows = []
            for i, j in slices:
                rows.append(np.mean(imputed_data[i:j], axis=0))
            return np.array(rows)

        true_row_values = data_without_nans.iloc[row_inds].values[:, 2:]

        diffs = true_row_values - impute_rows(imputed_data_methyLimp)
        add_metrics(diffs, rmse_lst, mae_lst)
        
        diffs = true_row_values - impute_rows(imputed_data_nbp.values[:, 2:].astype(np.float64))
        add_metrics(diffs, rmse_lst, mae_lst)

        diffs = true_row_values - impute_rows(imputed_data_people_mean.values[:, 2:].astype(np.float64))
        add_metrics(diffs, rmse_lst, mae_lst)

        diffs = true_row_values - impute_rows(imputed_data_cytosine_mean.values[:, 2:].astype(np.float64))
        add_metrics(diffs, rmse_lst, mae_lst)
    
    return rmse_lst, mae_lst


def combo_diffferent_missing_percent(filename, eps_methyLImp, eps_nbp1, eps_nbp2, p_arr=np.linspace(0.01, 0.2, 5)):
    """
    @param filename: path to the file with data
    @param eps_methyLImp: half neighborhood size for the methyLImp
    @param eps_nbp1: half neighborhood size for the nbp in the first part of the algorithm
    @param eps_nbp2: half neighborhood size for the nbp in the second part of the algorithm
    @param N_iter: number of iterations for each set of hyperparameters
    @param p: nans fraction
    @return: lists with rmse and mae for all sets of hyperparameters
    """
    data  = pd.read_csv(filename, sep='\t')
    mask = np.invert(np.any(np.isnan(data.values[:, 2:].astype(np.float64)), axis=1))
    data_without_nans = data[mask].reset_index(drop=True)
    
    n, _ = data_without_nans.shape

    row_inds = np.random.choice(np.arange(n), size=n // 1000, replace=False)
    data_without_rows = data_without_nans.drop(row_inds).reset_index(drop=True)
    true_row_values = data_without_nans.iloc[row_inds].values[:, 2:]

    n, m = data_without_rows.shape
    

    slices = []

    chromosome_coords = utils.get_chromosome_coords(data_without_rows.chromosome)
    positions = data_without_rows.position.values

    for row in row_inds:
        c = data_without_nans.chromosome[row]
        p = data_without_nans.position[row]
        i, j = chromosome_coords[c]
        l, r = utils.find_slice(p, eps_nbp2, positions[i:j])
        slices.append((l + i, r + i))

    rows_to_impute = np.zeros(len(data_without_rows), dtype=bool)
    for l, r in slices:
        rows_to_impute[l:r] = True

    def impute_rows(imputed_data):
            rows = []
            for i, j in slices:
                rows.append(np.mean(imputed_data[i:j], axis=0))
            return np.array(rows)
        
    rmse_methyLImp_lst = []
    mae_methyLImp_lst = []

    rmse_nbp_lst = []
    mae_nbp_lst = []

    for p in p_arr:
        print(f"missing values percent: {p}")
        
        nans_mask = np.random.choice([True, False], size=(n, m-2), p=[p, 1 - p])
        nans_mask = np.concatenate((np.zeros((n, 2), dtype=bool), nans_mask), axis=1)
        data_with_nans = data_without_rows.copy()
        data_with_nans[nans_mask] = np.nan

        impute_nans_mask = nans_mask & rows_to_impute[:, None]

        imputed_data_methyLImp = utils.methyLImp(data_with_nans, eps=eps_methyLImp, rows=rows_to_impute & np.any(nans_mask, axis=1))
        
        diffs_methyLImp = true_row_values - impute_rows(imputed_data_methyLImp)
        rmse_methyLImp_lst.append(np.mean(diffs_methyLImp ** 2) ** 0.5)
        mae_methyLImp_lst.append(np.mean(np.abs(diffs_methyLImp)))
        

        imputed_data_nbp = utils.impute_nbp(data_with_nans, eps=eps_nbp1, impute_positions=impute_nans_mask[:, 2:])

        diffs_nbp = true_row_values - impute_rows(imputed_data_nbp.values[:, 2:].astype(np.float64))
        rmse_nbp_lst.append(np.mean(diffs_nbp ** 2) ** 0.5)
        mae_nbp_lst.append(np.mean(np.abs(diffs_nbp)))

    return rmse_methyLImp_lst, mae_methyLImp_lst, rmse_nbp_lst, mae_nbp_lst
