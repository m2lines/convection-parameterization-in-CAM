"""Load a pytorch model and convert it to TorchScript."""
from typing import Optional
import torch

# FPTLIB-TODO
# Add a module import with your model here:
# This example assumes the model architecture is in an adjacent module `my_ml_model.py`
import models


def script_to_torchscript(
    model: torch.nn.Module, filename: Optional[str] = "scripted_model.pt"
) -> None:
    """
    Save pyTorch model to TorchScript using scripting.

    Parameters
    ----------
    model : torch.NN.Module
        a pyTorch model
    filename : str
        name of file to save to
    """
    # FIXME: torch.jit.optimize_for_inference() when PyTorch issue #81085 is resolved
    scripted_model = torch.jit.script(model)
    # print(scripted_model.code)
    scripted_model.save(filename)


def trace_to_torchscript(
    model: torch.nn.Module,
    dummy_input: torch.Tensor,
    filename: Optional[str] = "traced_model.pt",
) -> None:
    """
    Save pyTorch model to TorchScript using tracing.

    Parameters
    ----------
    model : torch.NN.Module
        a pyTorch model
    dummy_input : torch.Tensor
        appropriate size Tensor to act as input to model
    filename : str
        name of file to save to
    """
    # FIXME: torch.jit.optimize_for_inference() when PyTorch issue #81085 is resolved
    traced_model = torch.jit.trace(model, dummy_input)
    # traced_model.save(filename)
    frozen_model = torch.jit.freeze(traced_model)
    ## print(frozen_model.graph)
    ## print(frozen_model.code)
    frozen_model.save(filename)


def load_torchscript(filename: Optional[str] = "torchscript_model.pt") -> torch.nn.Module:
    """
    Load a TorchScript from file.

    Parameters
    ----------
    filename : str
        name of file containing TorchScript model
    """
    model = torch.jit.load(filename)

    return model


if __name__ == "__main__":
    # =====================================================
    # Load model and prepare for saving
    # =====================================================

    # FPTLIB-TODO
    # Load a pre-trained PyTorch model
    # Insert code here to load your model as `trained_model`.
    # This example assumes my_ml_model has a method `initialize` to load
    # architecture, weights, and place in inference mode

    # models.endow_with_netcdf_params(trained_model, weights_file)
    # weights_file = "qobsTTFFFFFTF30FFTFTF30TTFTFTFFF80FFTFTTF2699FFFF_X01_no_qp_no_adv_surf_F_Tin_qin_disteq_O_Trad_rest_Tadv_qadv_qout_qsed_RESCALED_7epochs_no_drop_REAL_NN_layers5in61out148_BN_F_te70.nc"
    # weights_file = "weights.pkl"
    # trained_model = models.ANN().initialize(weights_file, True)
    trained_model = models.ANN().load("nn_state.pt")
    # for name, layer in trained_model.named_children():
    #     print(layer.weight)
    #     break

    # Switch off specific layers/parts of the model that behave
    # differently during training and inference.
    # This may have been done by the user already, so just make sure here.
    trained_model.eval()

    # =====================================================
    # Prepare dummy input and check model runs
    # =====================================================

    # FPTLIB-TODO
    # Generate a dummy input Tensor `dummy_input` to the model of appropriate size.
    # This example assumes two inputs of size (512x40) and (512x1)
    trained_model_dummy_input = torch.ones((61), dtype=torch.float32)

    # FPTLIB-TODO
    # Uncomment the following lines to save for inference on GPU (rather than CPU):
    # device = torch.device('cuda')
    # trained_model = trained_model.to(device)
    # trained_model.eval()
    # trained_model_dummy_input = trained_model_dummy_input.to(device)

    # FPTLIB-TODO
    # Run model for dummy inputs
    # If something isn't working This will generate an error
    trained_model_dummy_output = trained_model(trained_model_dummy_input)

    # =====================================================
    # Save model
    # =====================================================

    # FPTLIB-TODO
    # Set the name of the file you want to save the torchscript model to:
    saved_ts_filename = "torchscript_model.pt"

    # FPTLIB-TODO
    # Save the pytorch model using either scripting (recommended where possible) or tracing
    # -----------
    # Scripting
    # -----------
    script_to_torchscript(trained_model, filename=saved_ts_filename)

    # -----------
    # Tracing
    # -----------
    # trace_to_torchscript(trained_model, trained_model_dummy_input, filename=saved_ts_filename)

    # =====================================================
    # Check model saved OK
    # =====================================================

    # Load torchscript and run model as a test
    # FPTLIB-TODO
    # Scale inputs as above and, if required, move inputs and mode to GPU
    trained_model_dummy_input = 2.0 * trained_model_dummy_input
    trained_model_testing_output = trained_model(trained_model_dummy_input)
    ts_model = load_torchscript(filename=saved_ts_filename)
    ts_model_output = ts_model(trained_model_dummy_input)

    if torch.all(ts_model_output.eq(trained_model_testing_output)):
        print("Saved TorchScript model working as expected in a basic test.")
        print("Users should perform further validation as appropriate.")
    else:
        raise RuntimeError(
            "Saved Torchscript model is not performing as expected.\n"
            "Consider using scripting if you used tracing, or investigate further."
        )
