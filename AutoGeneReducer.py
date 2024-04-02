import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from keras.models import Model, load_model
from keras.layers import Input, Dense, Dropout

def load_data(filepath):
    """
    Load data from a file using Pandas. This function automatically detects the file type
    based on its extension and chooses the appropriate Pandas function to load the data.
    Supported file types include CSV, XLS and  XLSX. For CSV and Excel files, the first column
    will be used as the index of the DataFrame if 'index_col=0' is applicable.

    Parameters:
    filepath (str): Path to the data file (i.e. Absolute or relative path to the data file).

    Returns:
    DataFrame: Loaded data as a Pandas DataFrame.

    Raises:
    - ValueError: If the file format is not supported.

    Example:
    >>> data = load_data('gene_expression.csv')
    This will load a CSV file containing gene expression data into a DataFrame, assuming the first
    column is an index column which containes the gene names/Ids.
    
    """
    _, file_extension = os.path.splitext(filepath)
    if file_extension.lower() in ['.csv']:
        data = pd.read_csv(filepath, index_col=0)
    elif file_extension.lower() in ['.xls', '.xlsx']:
        data = pd.read_excel(filepath, index_col=0)
    else:
        raise ValueError("Unsupported file type: Please provide a CSV, XLS or XLSX file.")
    print("Data loading complete.")
    return data

def build_ae(original_dim, params):
    """
    Constructs an autoencoder model using Keras based on the specified parameters. The model architecture is defined by the parameters
    passed. This includes the dimensionality of the input data, the size and number of layers, and the dropout rate to prevent overfitting.
    
    Parameters:
    original_dim (int): Dimensionality of the input data.
    - params (dict): Dictionary containing model hyperparameters, including:
        - first_layer_size (int): Size of the first encoding layer.
        - depth (int): Number of layers in the encoding part.
        - latent_dim (int): Size of the latent space.
        - dropout_rate (float): Dropout rate applied to each layer.
        - activation (str): Activation function used in the layers.
    
    Returns:
    Model: Compiled autoencoder model.
    """
    # Input and encoding layers
    input_layer = Input(shape=(original_dim,))
    encoded = Dense(params['first_layer_size'], activation=params['activation'])(input_layer)
    encoded = Dropout(params['dropout_rate'])(encoded)
    # Adding additional encoding layers based on depth
    for _ in range(params['depth'] - 1):
        encoded = Dense(params['first_layer_size'], activation=params['activation'])(encoded)
        encoded = Dropout(params['dropout_rate'])(encoded)
    # Latent space
    encoded = Dense(params['latent_dim'], activation=params['activation'])(encoded)
    # Decoding layer
    decoded = Dense(original_dim, activation=params['activation'])(encoded)
    # Model definition
    autoencoder = Model(input_layer, decoded)
    # Model compilation
    autoencoder.compile(optimizer='adam', loss='mse')
    return autoencoder
def train_ae(autoencoder, X_train, X_val, params, shuffle):
    """
    Trains the autoencoder model on the provided data.

    Parameters:
    autoencoder (Model): The autoencoder model to be trained.
    X_train (DataFrame): Training data.
    X_val (DataFrame): Validation data.
    params (dict): Dictionary containing training parameters.
    shuffle (bool): Whether to shuffle the training data before each epoch.

    Returns:
    Tuple(Model, History): The trained model and the training history object.
    """
    history = autoencoder.fit(X_train, X_train,
                              epochs=params['epochs'],
                              batch_size=params['batch_size'],
                              shuffle=shuffle,
                              validation_data=(X_val, X_val))
    return autoencoder, history

def plot_loss(history):
    """
    Plots the training and validation loss from the model's history.
    """
    plt.figure()
    plt.plot(history.history['loss'], label='Train Loss')
    plt.plot(history.history['val_loss'], label='Validation Loss')
    plt.title('Loss over Epochs')
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    plt.legend()
    plt.show()

def plot_reconstruction_error(autoencoder, X_val):
    """
    Plots the reconstruction error on the validation dataset.
    """
    reconstructions = autoencoder.predict(X_val)
    reconstruction_error = np.mean(np.abs(reconstructions - X_val), axis=1)
    plt.figure()
    plt.hist(reconstruction_error, bins=50)
    plt.title('Reconstruction Error Histogram')
    plt.xlabel('Reconstruction Error')
    plt.ylabel('Number of Samples')
    plt.show()

def main():
    """
    Main function to run the autoencoder pipeline.
    """
    parser = argparse.ArgumentParser(description="Run and train the autoencoder for dimensionality reduction and Feature extraction.")
    parser.add_argument('data_filepath', type=str, help='The file path of the data CSV,XLS or XLSX.')
    parser.add_argument('--val_data_filepath', type=str, default=None, help='Optional: The file path of the validation data CSV,XLS, XLSX. If not provided, a split from the training data will be used.')
    parser.add_argument('--model_path', type=str, default=None, help='Optional: The file path of the pre-trained model to load.')
    parser.add_argument('--latent_dim', type=int, default=20, help='Dimension of the latent space. Default is 20.')
    parser.add_argument('--depth', type=int, default=2, help='Depth of the neural network (number of layers). Default is 2.')
    parser.add_argument('--first_layer_size', type=int, default=500, help='Size of the first layer. Default is 500.')
    parser.add_argument('--dropout_rate', type=float, default=0.1, help='Dropout rate for training the neural network. Default is 0.1.')
    parser.add_argument('--epochs', type=int, default=200, help='Number of epochs to train the model. Default is 200.')
    parser.add_argument('--batch_size', type=int, default=80, help='Batch size for training the model. Default is 80.')
    parser.add_argument('--activation', type=str, default='relu', help="Activation function for the neural network layers. Default is 'relu'.")
    parser.add_argument('--test_size', type=float, default=0.25, help='Proportion of the dataset to include in the test split. Default is 0.25.')
    parser.add_argument('--shuffle', type=bool, default=True, help='Whether to shuffle the training data before each epoch. Default is True.')
    parser.add_argument('--plot_loss', action='store_true', help='Plot the loss curves after training.')
    parser.add_argument('--plot_reconstruction_error', action='store_true', help='Plot the reconstruction error histogram after training.')
    args = parser.parse_args()

    # Load training data
    data = load_data(args.data_filepath)
    original_dim = data.shape[0]  # Assuming features are in columns

    # Check if separate validation data is provided
    if args.val_data_filepath:
        X_val = load_data(args.val_data_filepath)
        X_train = data
    else:
        # Splitting into training and validation if no separate validation data provided
        data_transposed = data.transpose()
        X_train, X_val = train_test_split(data_transposed, test_size=args.test_size, random_state=42)

    params = {
        'latent_dim': args.latent_dim,
        'depth': args.depth,
        'first_layer_size': args.first_layer_size,
        'dropout_rate': args.dropout_rate,
        'epochs': args.epochs,
        'batch_size': args.batch_size,
        'activation': args.activation
    }

    # Model setup and training
    if args.model_path and os.path.exists(args.model_path):
        autoencoder = load_model(args.model_path)
        print("Loaded pre-trained model from:", args.model_path)
    else:
        autoencoder = build_ae(original_dim, params)
        print("Built new model.")

    trained_autoencoder, history = train_ae(autoencoder, X_train, X_val, params, args.shuffle)

    # Save the trained model
    trained_autoencoder.save('autoencoder_model.keras')
    print("Autoencoder model saved as 'autoencoder_model.keras'.")

    # Plotting
    if args.plot_loss:
        plot_loss(history)
    if args.plot_reconstruction_error:
        plot_reconstruction_error(trained_autoencoder, X_val)

if __name__ == "__main__":
    main()
