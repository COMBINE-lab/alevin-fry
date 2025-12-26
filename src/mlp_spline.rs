use anyhow::{Context, Result};
use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize};
// use std::error::Error;
use tch::nn::ModuleT;
use tch::no_grad;
use tch::{nn, Device,Kind, Tensor};
use std::fs::File;
use std::io::BufReader;
use serde_json::Value;

#[derive(Deserialize, Serialize)]
pub struct MLPParams {
    weights: Vec<Vec<Vec<f64>>>, // Deserialize as Vec<Vec<f64>> first
    biases: Vec<Vec<f64>>,       // Deserialize as Vec<f64> first
    activation: String,          // ReLU, Common activation function for all layers
    output_activation: String,   // logistic, Activation function for the output layer
}
pub struct MLPParamsND {
    weights: Vec<Array2<f64>>, // Converted to ndarray types
    biases: Vec<Array1<f64>>,  // Converted to ndarray types
    activation: String,
    output_activation: String,
}

impl MLPParams {
    fn to_ndarray(self) -> Result<MLPParamsND> {
        let weights = self
            .weights
            .into_iter()
            .map(|w| {
                Array2::from_shape_vec((w.len(), w[0].len()), w.into_iter().flatten().collect())
                    .with_context(|| "Failed to convert weights to ndarray")
            })
            .collect::<Result<Vec<Array2<f64>>>>()?;

        let biases = self
            .biases
            .into_iter()
            .map(|b| {
                Array1::from_shape_vec(b.len(), b)
                    .with_context(|| "Failed to convert biases to ndarray")
            })
            .collect::<Result<Vec<Array1<f64>>>>()?;

        Ok(MLPParamsND {
            weights,
            biases,
            activation: self.activation,
            output_activation: self.output_activation,
        })
    }
}

pub fn load_mlp_params(file_path: &str) -> Result<MLPParamsND> {
    let file = std::fs::File::open(file_path)
        .with_context(|| format!("Failed to open file at path: {}", file_path))?;
    let params: MLPParams = serde_json::from_reader(file)
        .with_context(|| "Failed to parse MLP parameters from JSON")?;
    params.to_ndarray()
}

pub fn create_mlp_with_tch(params: &MLPParamsND) -> Result<nn::Sequential> {
    let vs = nn::VarStore::new(Device::Cpu);
    let mut net = nn::seq();

    // Adjusted to reflect the transposed weights
    for (i, (weights, biases)) in params.weights.iter().zip(&params.biases).enumerate() {
        let output_size = weights.shape()[0]; // After transpose
        let input_size = weights.shape()[1];

        let mut layer = nn::linear(
            &vs.root(),
            input_size as i64,  // in_features
            output_size as i64, // out_features
            Default::default(),
        );

        no_grad(|| {
            // No need to transpose if weights are already transposed in Python
            layer.ws.copy_(
                &Tensor::from_slice(weights.as_slice().unwrap())
                    .reshape([output_size as i64, input_size as i64])
                    .to_kind(Kind::Float),
            );

            // Copy biases
            if let Some(ref mut bias) = layer.bs {
                bias.copy_(&Tensor::from_slice(biases.as_slice().unwrap()).to_kind(Kind::Float));
            }
        });

        net = net.add(layer);

        // Activation function
        if i < params.weights.len() - 1 {
            if params.activation.to_lowercase() == "relu" {
                net = net.add_fn(|x| x.relu());
            } else {
                return Err(anyhow::anyhow!(
                    "Unsupported activation: {}",
                    params.activation
                ));
            }
        }
    }

    // Output activation
    if params.output_activation.to_lowercase() == "logistic" {
        net = net.add_fn(|x| x.sigmoid()); // Use PyTorch's sigmoid
    } else {
        return Err(anyhow::anyhow!(
            "Unsupported output activation: {}",
            params.output_activation
        ));
    }

    Ok(net)
}

pub fn predict_with_tch(net: &nn::Sequential, input: Array2<f32>) -> Result<Array1<f64>> {
    let input_tensor = Tensor::from_slice(input.as_slice().unwrap())
        .reshape([input.shape()[0] as i64, input.shape()[1] as i64]);
    // Perform prediction
    let output_tensor = net.forward_t(&input_tensor, false);

    // Convert Tensor back to a 1D array of f64 (positive class scores)
    let output_vec: Vec<f64> = output_tensor
        .view(-1) // Flatten the tensor to 1D
        .shallow_clone() // Clone the tensor
        .try_into()
        .expect("Failed to convert Tensor to Vec<f64>");

    let output_array = Array1::from_vec(output_vec); // Create 1D Array1

    Ok(output_array)
}
// fn main() -> Result<()> {
//     // Load MLP parameters from JSON
//     let params = load_mlp_params("/fs/nexus-projects/sc_frag_len/nextflow/umi_level/my_clean_custom_alevin_fry/speed_up_PUG_forseti/convert_to_rust/mlp_params.json")?;

//     // Create MLP model using Tch
//     let mlp = create_mlp_with_tch(&params)?;

//     // Prepare input (one-hot encoded data with 150 columns)
//     let random_array: Array2<f32> = Array2::random((50, 150), Uniform::new(0.0, 1.0));

//     // Predict
//     let predictions = predict_with_tch(&mlp, random_array)?;

//     // Print predictions
//     println!("Predictions: {:?}", predictions);

//     Ok(())
// }


pub fn load_spline_lookup_table(file_path: &str) -> anyhow::Result<Array1<f64>> {
    // Open the file
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    // Parse the JSON
    let json_data: Value = serde_json::from_reader(reader)?;

    // Extract the "y" field as an array
    if let Some(y_values) = json_data["y"].as_array() {
        // Convert JSON array to Vec<f64>
        let y_vec: Vec<f64> = y_values
            .iter()
            .map(|v| v.as_f64().unwrap_or(0.0)) // Ensure conversion
            .collect();

        // Convert Vec<f64> to ndarray::Array1
        Ok(Array1::from(y_vec))
    } else {
        Err(anyhow::anyhow!("Missing 'y' field in JSON"))
    }
}