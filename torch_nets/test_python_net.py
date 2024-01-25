import models as nn_model

from torch import ones  # pylint: disable=no-name-in-module

if __name__ == "__main__":

    test_model = nn_model.ANN()

    nn_model.endow_with_netcdf_params(test_model, "qobsTTFFFFFTF30FFTFTF30TTFTFTFFF80FFTFTTF2699FFFF_X01_no_qp_no_adv_surf_F_Tin_qin_disteq_O_Trad_rest_Tadv_qadv_qout_qsed_RESCALED_7epochs_no_drop_REAL_NN_layers5in61out148_BN_F_te70.nc")

    test_model.eval()

    output = test_model(ones(61))

    print(output)

    with open("nn_out.txt", "w", encoding="utf_8") as file:
        for val in output:
            file.write(f"{val.item()}\n")
