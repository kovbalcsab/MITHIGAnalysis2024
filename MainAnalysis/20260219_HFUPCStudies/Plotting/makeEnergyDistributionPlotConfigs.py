import json
import os


pt_y_bins = [[2, 5, -2, 2]]

is_gamma_n_array = [1]

config_file_tag = "original_MC_HFEMaxPlus_forest_comparison"
plotting_dir = f"Plotting/EnergyDistributionPlots/{config_file_tag}"
plot_config_file_name = f"Plotting/energyDistributionPlotConfigs/PlotConfig_{config_file_tag}.json"

inputs = [
    #Unfolding nominal:
    #{"SampleDir": "TemplateFitting/output_unfolding/data_unfolding_ybinned", "FileName": "unfolding_HFEMaxPlus_forest_pt2-5_y-2-2.root", "SampleLabel": "-2-2", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/data_unfolding_ybinned", "FileName": "unfolding_HFEMaxPlus_forest_pt2-5_y-2--1.root", "SampleLabel": "-2--1", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/data_unfolding_ybinned", "FileName": "unfolding_HFEMaxPlus_forest_pt2-5_y-1-0.root", "SampleLabel": "-1-0", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/data_unfolding_ybinned", "FileName": "unfolding_HFEMaxPlus_forest_pt2-5_y0-1.root", "SampleLabel": "0-1", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/data_unfolding_ybinned", "FileName": "unfolding_HFEMaxPlus_forest_pt2-5_y1-2.root", "SampleLabel": "1-2", "AddKinematicRange": 0}
    #{"SampleDir": "TemplateFitting/output_unfolding/data_unfolding_ybinned", "FileName": "unfolding_HFEMaxMinus_forest_pt2-5_y-2-2.root", "SampleLabel": "-2-2", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/data_unfolding_ybinned", "FileName": "unfolding_HFEMaxMinus_forest_pt2-5_y-2--1.root", "SampleLabel": "-2--1", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/data_unfolding_ybinned", "FileName": "unfolding_HFEMaxMinus_forest_pt2-5_y-1-0.root", "SampleLabel": "-1-0", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/data_unfolding_ybinned", "FileName": "unfolding_HFEMaxMinus_forest_pt2-5_y0-1.root", "SampleLabel": "0-1", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/data_unfolding_ybinned", "FileName": "unfolding_HFEMaxMinus_forest_pt2-5_y1-2.root", "SampleLabel": "1-2", "AddKinematicRange": 0}
    #{"SampleDir": "TemplateFitting/output_unfolding/MC_unfolding_ybinned", "FileName": "unfoldingMC_HFEMaxPlus_forest_pt2-5_y-2-2.root", "SampleLabel": "-2-2", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/MC_unfolding_ybinned", "FileName": "unfoldingMC_HFEMaxPlus_forest_pt2-5_y-2--1.root", "SampleLabel": "-2--1", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/MC_unfolding_ybinned", "FileName": "unfoldingMC_HFEMaxPlus_forest_pt2-5_y-1-0.root", "SampleLabel": "-1-0", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/MC_unfolding_ybinned", "FileName": "unfoldingMC_HFEMaxPlus_forest_pt2-5_y0-1.root", "SampleLabel": "0-1", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/MC_unfolding_ybinned", "FileName": "unfoldingMC_HFEMaxPlus_forest_pt2-5_y1-2.root", "SampleLabel": "1-2", "AddKinematicRange": 0}

    #SVD:
    #{"SampleDir": "TemplateFitting/output_unfolding/SVD_SVD_MC_unfolding_ybinned", "FileName": "unfoldingMC_HFEMaxPlus_forest_SVD_pt2-5_y-2-2.root", "SampleLabel": "-2-2", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/SVD_SVD_MC_unfolding_ybinned", "FileName": "unfoldingMC_HFEMaxPlus_forest_SVD_pt2-5_y-2--1.root", "SampleLabel": "-2--1", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/SVD_SVD_MC_unfolding_ybinned", "FileName": "unfoldingMC_HFEMaxPlus_forest_SVD_pt2-5_y-1-0.root", "SampleLabel": "-1-0", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/SVD_SVD_MC_unfolding_ybinned", "FileName": "unfoldingMC_HFEMaxPlus_forest_SVD_pt2-5_y0-1.root", "SampleLabel": "0-1", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/SVD_SVD_MC_unfolding_ybinned", "FileName": "unfoldingMC_HFEMaxPlus_forest_SVD_pt2-5_y1-2.root", "SampleLabel": "1-2", "AddKinematicRange": 0}

    #FFT:
    #{"SampleDir": "TemplateFitting/output_unfolding/FFT_FFT_MC_deconvolution_ybinned", "FileName": "deconvolveFFTMC_HFEMaxPlus_forest_pt2-5_y-2-2.root", "SampleLabel": "-2-2", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/FFT_FFT_MC_deconvolution_ybinned", "FileName": "deconvolveFFTMC_HFEMaxPlus_forest_pt2-5_y-2--1.root", "SampleLabel": "-2--1", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/FFT_FFT_MC_deconvolution_ybinned", "FileName": "deconvolveFFTMC_HFEMaxPlus_forest_pt2-5_y-1-0.root", "SampleLabel": "-1-0", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/FFT_FFT_MC_deconvolution_ybinned", "FileName": "deconvolveFFTMC_HFEMaxPlus_forest_pt2-5_y0-1.root", "SampleLabel": "0-1", "AddKinematicRange": 0},
    #{"SampleDir": "TemplateFitting/output_unfolding/FFT_FFT_MC_deconvolution_ybinned", "FileName": "deconvolveFFTMC_HFEMaxPlus_forest_pt2-5_y1-2.root", "SampleLabel": "1-2", "AddKinematicRange": 0}

    #Measured:
    {"SampleDir": "TemplateFitting/output_unfolding/FFT_FFT_MC_deconvolution_ybinned", "FileName": "deconvolveFFTMC_HFEMaxPlus_forest_pt2-5_y-2-2.root", "SampleLabel": "-2-2", "AddKinematicRange": 0},
    {"SampleDir": "TemplateFitting/output_unfolding/FFT_FFT_MC_deconvolution_ybinned", "FileName": "deconvolveFFTMC_HFEMaxPlus_forest_pt2-5_y-2--1.root", "SampleLabel": "-2--1", "AddKinematicRange": 0},
    {"SampleDir": "TemplateFitting/output_unfolding/FFT_FFT_MC_deconvolution_ybinned", "FileName": "deconvolveFFTMC_HFEMaxPlus_forest_pt2-5_y-1-0.root", "SampleLabel": "-1-0", "AddKinematicRange": 0},
    {"SampleDir": "TemplateFitting/output_unfolding/FFT_FFT_MC_deconvolution_ybinned", "FileName": "deconvolveFFTMC_HFEMaxPlus_forest_pt2-5_y0-1.root", "SampleLabel": "0-1", "AddKinematicRange": 0},
    {"SampleDir": "TemplateFitting/output_unfolding/FFT_FFT_MC_deconvolution_ybinned", "FileName": "deconvolveFFTMC_HFEMaxPlus_forest_pt2-5_y1-2.root", "SampleLabel": "1-2", "AddKinematicRange": 0}
]

energy_distributions = [
    #Unfolding nominal:
    #{"VarName": "hUnfolded_iter7", "PlotBaseName": "HFEMaxPlus_data_unfolding_iter7_energy_comparison_zoom.pdf",
    # "XTitle": "HF~E_{max}^{+}~[GeV]", "YTitle": "Normalized~Entries", "PlotTitle": "HF~E_{max}^{+}~distribution~comparison",
    # "LogY": 1, "XLims": [0, 25], "YLims": [1e-4, 1.2], "Normalize": 1, "RebinFactor": 4, "IsTreeVar": 0, "DoRatio": 1},
    #{"VarName": "hUnfolded_iter7", "PlotBaseName": "HFEMaxPlus_data_unfolding_iter7_energy_comparison_full.pdf",
    # "XTitle": "HF~E_{max}^{+}~[GeV]", "YTitle": "Normalized~Entries", "PlotTitle": "HF~E_{max}^{+}~distribution~comparison",
    # "LogY": 1, "XLims": [0, 200], "YLims": [1e-4, 1.2], "Normalize": 1, "RebinFactor": 4, "IsTreeVar": 0, "DoRatio": 1}
    #{"VarName": "hUnfolded_iter7", "PlotBaseName": "HFEMaxMinus_data_unfolding_iter7_energy_comparison_zoom.pdf",
    # "XTitle": "HF~E_{max}^{-}~[GeV]", "YTitle": "Normalized~Entries", "PlotTitle": "HF~E_{max}^{-}~distribution~comparison",
    # "LogY": 1, "XLims": [0, 25], "YLims": [1e-4, 1.2], "Normalize": 1, "RebinFactor": 4, "IsTreeVar": 0, "DoRatio": 1},
    #{"VarName": "hUnfolded_iter7", "PlotBaseName": "HFEMaxMinus_data_unfolding_iter7_energy_comparison_full.pdf",
    # "XTitle": "HF~E_{max}^{-}~[GeV]", "YTitle": "Normalized~Entries", "PlotTitle": "HF~E_{max}^{-}~distribution~comparison",
    # "LogY": 1, "XLims": [0, 200], "YLims": [1e-4, 1.2], "Normalize": 1, "RebinFactor": 4, "IsTreeVar": 0, "DoRatio": 1}

    #SVD:
    #{"VarName": "hUnfolded_iter6", "PlotBaseName": "HFEMaxPlus_data_unfolding_iter6_energy_comparison_zoom.pdf",
    # "XTitle": "HF~E_{max}^{+}~[GeV]", "YTitle": "Normalized~Entries", "PlotTitle": "HF~E_{max}^{+}~distribution~comparison",
    # "LogY": 1, "XLims": [0, 25], "YLims": [1e-4, 1.2], "Normalize": 1, "RebinFactor": 4, "IsTreeVar": 0, "DoRatio": 1},
    #{"VarName": "hUnfolded_iter6", "PlotBaseName": "HFEMaxPlus_data_unfolding_iter6_energy_comparison_full.pdf",
    # "XTitle": "HF~E_{max}^{+}~[GeV]", "YTitle": "Normalized~Entries", "PlotTitle": "HF~E_{max}^{+}~distribution~comparison",
    # "LogY": 1, "XLims": [0, 200], "YLims": [1e-4, 1.2], "Normalize": 1, "RebinFactor": 4, "IsTreeVar": 0, "DoRatio": 1}

    #FFT:
    #{"VarName": "hDeconvolved_iter4", "PlotBaseName": "HFEMaxPlus_data_unfolding_iter4_energy_comparison_zoom.pdf",
    # "XTitle": "HF~E_{max}^{+}~[GeV]", "YTitle": "Normalized~Entries", "PlotTitle": "HF~E_{max}^{+}~distribution~comparison",
    # "LogY": 1, "XLims": [0, 25], "YLims": [1e-4, 1.2], "Normalize": 1, "RebinFactor": 4, "IsTreeVar": 0, "DoRatio": 1},
    #{"VarName": "hDeconvolved_iter4", "PlotBaseName": "HFEMaxPlus_data_unfolding_iter4_energy_comparison_full.pdf",
    # "XTitle": "HF~E_{max}^{+}~[GeV]", "YTitle": "Normalized~Entries", "PlotTitle": "HF~E_{max}^{+}~distribution~comparison",
    # "LogY": 1, "XLims": [0, 200], "YLims": [1e-4, 1.2], "Normalize": 1, "RebinFactor": 4, "IsTreeVar": 0, "DoRatio": 1}
    
    #Measured:
    {"VarName": "hMeasured", "PlotBaseName": "HFEMaxPlus_dataMeasured_energy_comparison_zoom.pdf",
     "XTitle": "HF~E_{max}^{+}~[GeV]", "YTitle": "Normalized~Entries", "PlotTitle": "HF~E_{max}^{+}~distribution~comparison",
     "LogY": 1, "XLims": [0, 25], "YLims": [1e-4, 1.2], "Normalize": 1, "RebinFactor": 4, "IsTreeVar": 0, "DoRatio": 1},
    {"VarName": "hMeasured", "PlotBaseName": "HFEMaxPlus_dataMeasured_energy_comparison_full.pdf",
     "XTitle": "HF~E_{max}^{+}~[GeV]", "YTitle": "Normalized~Entries", "PlotTitle": "HF~E_{max}^{+}~distribution~comparison",
     "LogY": 1, "XLims": [0, 200], "YLims": [1e-4, 1.2], "Normalize": 1, "RebinFactor": 4, "IsTreeVar": 0, "DoRatio": 1}
]

plots = []
for is_gamma_n in is_gamma_n_array:
    for pt_min, pt_max, y_min, y_max in pt_y_bins:
        for dist in energy_distributions:
            curr_input_files = []
            curr_labels = []
            curr_is_tree = []
            curr_n_bins = []
            curr_var_names = []
            for item in inputs:
                curr_labels.append(item["SampleLabel"])
                if item["AddKinematicRange"]:
                    curr_input_files.append(
                        f"{item['SampleDir']}/pt{pt_min}-{pt_max}_y{y_min}-{y_max}_IsGammaN{is_gamma_n}/{item['FileName']}"
                    )
                else:
                    curr_input_files.append(f"{item['SampleDir']}/{item['FileName']}")
                curr_is_tree.append(str(dist["IsTreeVar"]))
                if dist["IsTreeVar"]:
                    curr_n_bins.append(str(dist["NBins"]))
                curr_var_names.append(dist["VarName"])

            plots.append({
                "PlotBaseName": dist["PlotBaseName"],
                "StatMatrixBaseName": dist["PlotBaseName"].replace(".pdf", "_Chi2PValueMatrix.pdf"),
                "WassersteinMatrixBaseName": dist["PlotBaseName"].replace(".pdf", "_WassersteinMatrix.pdf"),
                "JSDistanceMatrixBaseName": dist["PlotBaseName"].replace(".pdf", "_JSDistanceMatrix.pdf"),
                "HellingerMatrixBaseName": dist["PlotBaseName"].replace(".pdf", "_HellingerMatrix.pdf"),
                "PoissonDevianceMatrixBaseName": dist["PlotBaseName"].replace(".pdf", "_PoissonDevianceMatrix.pdf"),
                "InputFileNames": ",".join(curr_input_files),
                "VarNames": ",".join(curr_var_names),
                "IsTreeVar": ",".join(curr_is_tree),
                "Labels": ",".join(curr_labels),
                "xTitle": dist["XTitle"],
                "yTitle": dist["YTitle"],
                "plotTitle": dist["PlotTitle"],
                "logY": dist["LogY"],
                "xMin": dist["XLims"][0],
                "xMax": dist["XLims"][1],
                "yMin": dist["YLims"][0],
                "yMax": dist["YLims"][1],
                "normalizeToUnity": dist["Normalize"],
                "rebinFactor": dist["RebinFactor"],
                "StatRebin": dist.get("StatRebin", 1),
                "StatEpsilon": dist.get("StatEpsilon", 1e-12),
                "PoissonFloorMu": dist.get("PoissonFloorMu", 1e-9),
                "NBins": ",".join(curr_n_bins), 
                "doRatio": dist["DoRatio"],
                "MinDzeroPT": pt_min,
                "MaxDzeroPT": pt_max,
                #"MinDzeroY": y_min,
                #"MaxDzeroY": y_max,
                "IsGammaN": is_gamma_n
            })

data = {
    "plottingDir": plotting_dir,
    "Plots": plots
}

os.makedirs(os.path.dirname(plot_config_file_name), exist_ok=True)
with open(plot_config_file_name, "w") as out:
    json.dump(data, out, indent=4)

print(f"Wrote config: {plot_config_file_name}")
