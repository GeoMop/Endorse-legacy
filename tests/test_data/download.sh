#!/bin/bash

# download big test data
curl https://flow.nti.tul.cz/endorse_large_data/large_model.tar.gz --output large_model.tar.gz
tar -xvf large_model.tar.gz

curl https://flow.nti.tul.cz/endorse_large_data/large_model_local.tar.gz --output large_model_local.tar.gz
tar -xvf large_model_local.tar.gz
