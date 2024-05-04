import pickle

import pandas as pd


def predict_operon(model_file: str = "../data/model.pkl", data_for_predict: str = pd.DataFrame) -> pd.Series:
    """
    :param model_file: path to model
    :return:pd.Series with predicted operons
    """

    with open(model_file, "rb") as f:
        model = pickle.load(f)

    predictions = model.predict(data_for_predict)
    return pd.Series(predictions)
