import pandas as pd
from random import choice

df = pd.read_table("Elenco-comuni-italiani.csv", sep=",")

numeri_rimanenti = list(range(len(df)))
punteggio = 0
punteggio_max = 0

df[df.SIGLA == "".upper()]
df[df.provincia == "".capitalize()]

while True:
    casuale = choice(numeri_rimanenti)
    sigla = df.iloc[casuale].SIGLA
    if type(sigla) is float:
        sigla = "NA"
    comune = df.iloc[casuale].COMUNE
    ripartizione_geografica = df["Ripartizione geografica"].iloc[casuale]
    provincia = df.iloc[casuale].provincia
    regione = df.iloc[casuale].regione

    risposta = input(comune + " ")
    risposta = risposta.strip()

    sigla_risposta = df[df.SIGLA == risposta.upper()]
    provincia_risposta = df[df.provincia == risposta.capitalize()]

    if len(provincia_risposta) != 0:
        regione_risposta = str(provincia_risposta.regione.iloc[0])
        ripartizione_geografica_risposta = str(
            provincia_risposta["Ripartizione geografica"].iloc[0]
        )
    elif len(sigla_risposta) != 0:
        regione_risposta = str(sigla_risposta.regione.iloc[0])
        ripartizione_geografica_risposta = str(
            sigla_risposta["Ripartizione geografica"].iloc[0]
        )
    else:
        regione_risposta = None
        ripartizione_geografica_risposta = None

    if risposta.upper() == sigla or risposta.capitalize() == provincia:
        punteggio += 5
        print("Giusto!")
    elif regione_risposta == regione:
        punteggio += 3
        print(f"No, era {sigla}")
        print("Regione giusta")
    elif ripartizione_geografica_risposta == ripartizione_geografica:
        punteggio += 1
        print(f"No, era {sigla}")
        print("Ripartizione geografica giusta")
    else:
        print(f"No, era {sigla}")
    punteggio_max += 5
    print(f"punteggio: {punteggio}/{punteggio_max}")
    numeri_rimanenti.remove(casuale)

# print(sigla,comune,regione,provincia,ripartizione_geografica)
