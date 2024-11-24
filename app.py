import streamlit as st
import pandas as pd
import plotly.express as px
from Bio import SeqIO
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
        "AUG":"M", "UGG":"W", "UUU":"F", "UUC":"F",
        "UUA":"L", "UUG":"L", "CUU":"L", "CUC":"L",
        "CUA":"L", "CUG":"L", "AUU":"I", "AUC":"I",
        "AUA":"I", "GUU":"V", "GUC":"V", "GUA":"V",
        "GUG":"V", "UCU":"S", "UCC":"S", "UCA":"S",
        "UCG":"S", "CCU":"P", "CCC":"P", "CCA":"P",
        "CCG":"P", "ACU":"T", "ACC":"T", "ACA":"T",
        "ACG":"T", "GCU":"A", "GCC":"A", "GCA":"A",
        "GCG":"A", "UAU":"Y", "UAC":"Y", "CAU":"H",
        "CAC":"H", "CAA":"Q", "CAG":"Q", "AAU":"N",
        "AAC":"N", "AAA":"K", "AAG":"K", "GAU":"D",
        "GAC":"D", "GAA":"E", "GAG":"E", "UGU":"C",
        "UGC":"C", "UGA":"_", "UAA":"_", "UAG":"_",
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
input_sequence = st.sidebar.text_area(
    "Secuencia de ADN/ARN o proteína:",
    placeholder="Ejemplo: ATCGTTAGC o MVLTI...",
    height=150
)

# Botón para cargar ejemplos
if st.sidebar.button("Cargar Ejemplo"):
    input_sequence = "ATCGTTAGC"  # Ejemplo predefinido para ADN

# Subida de archivo FASTA
st.sidebar.subheader("Sube un archivo FASTA")
fasta_file = st.sidebar.file_uploader("Selecciona un archivo FASTA", type=["fasta"])

if fasta_file is not None:
    try:
        # Asegúrate de abrir el archivo en modo binario para evitar el error StreamModeError
        fasta_sequences = list(SeqIO.parse(fasta_file, "fasta"))
        
        if fasta_sequences:
            # Si se encuentra al menos una secuencia en el archivo FASTA
            input_sequence = str(fasta_sequences[0].seq)
            st.write("**Secuencia cargada desde el archivo FASTA:**")
            st.code(input_sequence)
        else:
            st.error("No se encontró ninguna secuencia en el archivo FASTA.")
    except Exception as e:
        st.error(f"Ocurrió un error al procesar el archivo FASTA: {e}")

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
        results = results.append({"Propiedad": "Contenido GC (%)", "Valor": gc_content}, ignore_index=True)

    if seq_type == "Proteína":
        results = results.append(
            {"Propiedad": "Residuos hidrofóbicos", "Valor": hydrophobic},
            ignore_index=True
        )
        results = results.append(
            {"Propiedad": "Residuos hidrofílicos", "Valor": hydrophilic},
            ignore_index=True
        )

    st.download_button(
        label="Descargar resultados como CSV",
        data=results.to_csv(index=False),
        file_name="resultados_bioinformaticos.csv",
        mime="text/csv"
    )
