import streamlit as st
import pandas as pd
import plotly.express as px
from Bio import SeqIO
import random
from collections import Counter

# Configuración inicial del Dashboard
st.set_page_config(page_title="Análisis Bioinformático", layout="centered")
st.title("Dashboard Bioinformático Interactivo")
st.sidebar.header("Opciones")

# Funciones auxiliares
def calcular_gc(sequence):
    """Calcula el contenido de GC en una secuencia de ADN."""
    gc_content = (sequence.count("G") + sequence.count("C")) / len(sequence) * 100
    return round(gc_content, 2)

def transcribir_adn(adn):
    """Transcribe una secuencia de ADN a ARN."""
    return adn.replace("T", "U")

def traducir_arn(arn):
    """Traduce una secuencia de ARN a proteína."""
    codigo_genetico = {
        "AUG": "M", "UGG": "W", "UUU": "F", "UUC": "F",
        "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L",
        "CUA": "L", "CUG": "L", "AUU": "I", "AUC": "I",
        "AUA": "I", "GUU": "V", "GUC": "V", "GUA": "V",
        "GUG": "V", "UCU": "S", "UCC": "S", "UCA": "S",
        "UCG": "S", "CCU": "P", "CCC": "P", "CCA": "P",
        "CCG": "P", "ACU": "T", "ACC": "T", "ACA": "T",
        "ACG": "T", "GCU": "A", "GCC": "A", "GCA": "A",
        "GCG": "A", "UAU": "Y", "UAC": "Y", "CAU": "H",
        "CAC": "H", "CAA": "Q", "CAG": "Q", "AAU": "N",
        "AAC": "N", "AAA": "K", "AAG": "K", "GAU": "D",
        "GAC": "D", "GAA": "E", "GAG": "E", "UGU": "C",
        "UGC": "C", "UGA": "_", "UAA": "_", "UAG": "_",
    }
    proteina = ""
    for i in range(0, len(arn) - 2, 3):
        codon = arn[i:i+3]
        proteina += codigo_genetico.get(codon, "X")
    return proteina

def visualizar_secuencia(sequence):
    """Visualiza una secuencia agrupada en bloques de 10 caracteres con colores."""
    bloques = [sequence[i:i+10] for i in range(0, len(sequence), 10)]
    colores = ["#FFDDC1", "#FFABAB", "#FFC3A0", "#D5AAFF"]
    html = "<div>"
    for i, bloque in enumerate(bloques):
        color = colores[i % len(colores)]
        html += f'<span style="background-color:{color};padding:4px;margin:2px;display:inline-block;">{bloque}</span>'
    html += "</div>"
    return html

def graficar_frecuencia(sequence):
    """Grafica la frecuencia de bases o residuos en la secuencia."""
    conteo = Counter(sequence)
    fig = px.bar(
        x=list(conteo.keys()),
        y=list(conteo.values()),
        title="Frecuencia de bases/residuos",
        labels={"x": "Base/Residuo", "y": "Frecuencia"},
        color=list(conteo.keys())
    )
    return fig

# Input de la secuencia
st.sidebar.subheader("Introduce tu secuencia")
input_type = st.sidebar.selectbox(
    "Selecciona el tipo de entrada:",
    ["Secuencia de ADN/ARN", "Secuencia de proteína", "ID o nombre de secuencia (NCBI)"]
)

if input_type == "ID o nombre de secuencia (NCBI)":
    input_sequence = st.sidebar.text_input(
        "Introduce un nombre de gen o ID (Ej. BRCA1):",
        placeholder="Ejemplo: NM_007294"
    )
    st.sidebar.write("Puedes consultar más detalles en [NCBI](https://www.ncbi.nlm.nih.gov/).")
else:
    input_sequence = st.sidebar.text_area(
        "Introduce tu secuencia:",
        placeholder="Ejemplo: ATCGTTAGC o MVLTI...",
        height=150
    )

# Botón para cargar ejemplos
if st.sidebar.button("Cargar Ejemplo"):
    if input_type == "Secuencia de ADN/ARN":
        input_sequence = "ATCGTTAGC"  # Ejemplo predefinido para ADN
    elif input_type == "Secuencia de proteína":
        input_sequence = "MVLTIHP"  # Ejemplo predefinido para proteína

# Validación y limpieza de la secuencia
if input_sequence:
    sequence = ''.join(filter(str.isalpha, input_sequence)).upper()
    st.write("**Secuencia procesada:**")
    st.code(sequence)

    # Identificar tipo de secuencia
    if set(sequence).issubset({"A", "T", "C", "G"}):
        seq_type = "ADN"
    elif set(sequence).issubset({"A", "U", "C", "G"}):
        seq_type = "ARN"
    else:
        seq_type = "Proteína"
    st.write(f"Tipo de secuencia detectado: **{seq_type}**")

    # Análisis
    length = len(sequence)
    st.write(f"Longitud de la secuencia: **{length}** bases o residuos")

    if seq_type == "ADN":
        gc_content = calcular_gc(sequence)
        st.write(f"Contenido GC: **{gc_content}%**")
        arn = transcribir_adn(sequence)
        st.write("**Transcripción a ARN:**")
        st.code(arn)

    if seq_type == "ARN":
        proteina = traducir_arn(sequence)
        st.write("**Traducción a Proteína:**")
        st.code(proteina)

    if seq_type == "Proteína":
        hydrophobic = sum(sequence.count(aa) for aa in "AILMFWYV")
        hydrophilic = sum(sequence.count(aa) for aa in "RKDENQ")

        # Gráfico de composición
        fig = px.pie(
            values=[hydrophobic, hydrophilic, length - hydrophobic - hydrophilic],
            names=["Hidrofóbicos", "Hidrofílicos", "Otros"],
            title="Composición de la proteína"
        )
        st.plotly_chart(fig)

    # Gráfico de frecuencia
    freq_fig = graficar_frecuencia(sequence)
    st.plotly_chart(freq_fig)

    # Visualización de la secuencia
    st.write("**Secuencia visualizada:**")
    st.markdown(visualizar_secuencia(sequence), unsafe_allow_html=True)

    # Exportar resultados
    results = pd.DataFrame({
        "Propiedad": ["Tipo de secuencia", "Longitud"],
        "Valor": [seq_type, length]
    })

    if seq_type == "ADN":
        new_row = pd.DataFrame([{"Propiedad": "Contenido GC (%)", "Valor": gc_content}])
        results = pd.concat([results, new_row], ignore_index=True)

    if seq_type == "Proteína":
        hydrophobic_row = pd.DataFrame([{"Propiedad": "Residuos hidrofóbicos", "Valor": hydrophobic}])
        hydrophilic_row = pd.DataFrame([{"Propiedad": "Residuos hidrofílicos", "Valor": hydrophilic}])
        results = pd.concat([results, hydrophobic_row, hydrophilic_row], ignore_index=True)

    st.download_button(
        label="Descargar resultados como CSV",
        data=results.to_csv(index=False),
        file_name="resultados_bioinformaticos.csv",
        mime="text/csv"
    )
