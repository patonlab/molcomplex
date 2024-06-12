import os

#ENV = "DEV"
ENV = "PROD"


## server
host = "0.0.0.0"
port = int(os.environ.get("PORT", 5000))


## info
app_name = "Molecular Complexity"
contacts = "https://ccas.nd.edu/"
code = "https://github.com/bobbypaton/molcomplex"
#tutorial = "https://towardsdatascience.com/web-development-with-python-dash-complete-tutorial-6716186e09b3"
fontawesome = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"

about = "Load your molecule to break its bonds"

## fs
#root = os.path.dirname(os.path.dirname(__file__)) + "/"