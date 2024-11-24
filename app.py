import streamlit as st
from Bio import Entrez, SeqIO
import pandas as pd
import plotly.express as px
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

# Función para obtener secuencia de NCBI
def obtener_secuencia_nucleotide(query, retmax=1):
    """Obtiene secuencias desde NCBI por nombre o ID de gen."""
    Entrez.email = "your-email@example.com"  # Cambia esto por tu email
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=retmax)
    record = Entrez.read(handle)
    id_gen = record["IdList"][0]  # Toma el primer gen de la búsqueda
    handle = Entrez.efetch(db="nucleotide", id=id_gen, rettype="gb", retmode="text")
    seq_record = SeqIO.read(handle, "genbank")
    return str(seq_record.seq)

# Input de la secuencia
st.sidebar.subheader("Buscar gen en NCBI")
gene = st.sidebar.text_input("Introduce un nombre de gen o ID (Ej. BRCA1):", "BRCA1")

# Obtener secuencia desde NCBI
if gene:
    try:
        st.sidebar.text("Buscando secuencia en NCBI...")
        sequence = obtener_secuencia_nucleotide(gene)
        st.write(f"**Secuencia obtenida para el gen {gene}:**")
        st.code(sequence)

        # Procesar la secuencia
        if sequence:
            seq_type = "ADN" if set(sequence).issubset({"A", "T", "C", "G"}) else "Desconocido"
            st.write(f"Tipo de secuencia detectado: **{seq_type}**")
            length = len(sequence)
            st.write(f"Longitud de la secuencia: **{length}** bases")

            if seq_type == "ADN":
                gc_content = calcular_gc(sequence)
                st.write(f"Contenido GC: **{gc_content}%**")
                arn = transcribir_adn(sequence)
                st.write("**Transcripción a ARN:**")
                st.code(arn)

            # Gráfico de frecuencia
            freq_fig = graficar_frecuencia(sequence)
            st.plotly_chart(freq_fig)

            # Exportar resultados
            results = pd.DataFrame({
                "Propiedad": ["Tipo de secuencia", "Longitud", "GC (%)"],
                "Valor": [seq_type, length, gc_content]
            })
            st.download_button(
                label="Descargar resultados como CSV",
                data=results.to_csv(index=False),
                file_name=f"resultados_{gene}.csv",
                mime="text/csv"
            )
    except Exception as e:
        st.error(f"Error al obtener la secuencia: {e}")
