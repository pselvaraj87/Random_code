import numpy as np
import json
import matplotlib.pyplot as plt
import pandas as pd


csvfile = 'C:\Users\pselvaraj\Desktop\Jaline_script\Lino_script\smeardata.csv'
smearing_values = pd.read_csv(csvfile, names=['log10_density_per_ml', 'log10_sigma'])

poly = np.polyfit(smearing_values.log10_density_per_ml,
                      smearing_values.log10_sigma,
                      deg=4)


def smear_density_values(true_densities):

    log10_true_densities = np.log10(true_densities)
    log10_sigma = np.polyval(poly, log10_true_densities + 3)
    log10_measured_densities = log10_true_densities + log10_sigma * np.random.normal(size=true_densities.size) / 2.  # 95% BCI
    measured_densities = np.power(10, log10_measured_densities)
    measured_densities[measured_densities < 1e-2] = 0

    measured_densities[np.isnan(measured_densities)] = 0
    measured_densities[true_densities < 2e-2] = 0
    return measured_densities


def smeared_density(positive_fields):

    uL_per_field = 0.5/200.0
    PfPR_smeared = -(1.0 / uL_per_field) * np.log(1 - positive_fields / 200.0)
    return PfPR_smeared


if __name__ == '__main__':
    survey_file = "MalariaSurveyJSONAnalyzer_Day_15_0.json"
    summary_file = "MalariaSummaryReport_Monthly_Report.json"
    with open(survey_file) as data_file:
        survey_data = json.load(data_file)
    with open(summary_file) as data_file:
        summary_data = json.load(data_file)

    true_asexual_density = np.array(survey_data['patient_array'][0]['true_asexual_parasites'])
    true_asexual_density_with_uncertainty = np.array(survey_data['patient_array'][0]['smeared_true_asexual_parasites'])
    true_gameto_density = np.array(survey_data['patient_array'][0]['true_gametocytes'])
    true_gameto_density_with_uncertainty = np.array(survey_data['patient_array'][0]['smeared_true_gametocytes'])

    plt.figure()
    plt.title('True parasite density')
    plt.plot(true_asexual_density, label='True Density')
    plt.plot(true_asexual_density_with_uncertainty, label='True Density with Uncertainty')
    plt.legend()
    plt.savefig('true_parasite_density.png')

    plt.figure()
    plt.title('True gametocyte density')
    plt.plot(true_gameto_density, label='True Density')
    plt.plot(true_gameto_density_with_uncertainty, label='True Density with Uncertainty')
    plt.legend()
    plt.savefig('true_gametocyte_density.png')

    plt.figure()
    plt.title('True vs. Corrected Asexual Parasite density')
    plt.scatter(true_asexual_density, true_asexual_density_with_uncertainty)
    plt.xlabel('True density')
    plt.ylabel('Corrected density')
    plt.xlim([0, 10000])
    plt.ylim([0, 10000])
    plt.savefig('scatter_parasite.png')

    plt.figure()
    plt.title('True vs. Corrected Gametocyte density')
    plt.scatter(true_gameto_density, true_gameto_density_with_uncertainty)
    plt.xlabel('True density')
    plt.ylabel('Corrected density')
    plt.ylim([0, 10000])
    plt.xlim([0, 10000])
    plt.savefig('scatter_gametocyte.png')