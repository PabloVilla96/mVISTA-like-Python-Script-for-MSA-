###### Este es el script que use para todas las figuras
###### Fue hecho con la herramienta de inteligenica artifcial Claude de Anthropic
###### Por Juan Pablo Villa Núñez, el 12 de noviembre de 2026


import os
import itertools
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from matplotlib.ticker import MultipleLocator

# Crear el directorio de salida si no existe
output_dir = r"D:\Maestría_DIC2025\Figuras_Tesis"
os.makedirs(output_dir, exist_ok=True)

# ==========================================
# 1. DICCIONARIO CENTRAL DE COORDENADAS Cgi
# ==========================================
c_prom, c_utr5, c_exon, c_intron, c_utr3 = '#ff9999', '#99ff99', '#80bfff', '#d9d9d9', '#df80ff'

estructuras_cgi = {
    'CgiPLT1/2a': [
        ('Promotor', 2000, c_prom), ("5' UTR", 90, c_utr5),
        ('Exón', 149, c_exon), ('Intrón', 267, c_intron),
        ('Exón', 441, c_exon), ('Intrón', 102, c_intron),
        ('Exón', 84, c_exon), ('Intrón', 147, c_intron),
        ('Exón', 10, c_exon), ('Intrón', 110, c_intron),
        ('Exón', 90, c_exon), ('Intrón', 117, c_intron),
        ('Exón', 72, c_exon), ('Intrón', 85, c_intron),
        ('Exón', 130, c_exon), ('Intrón', 500, c_intron),
        ('Exón', 785, c_exon), ("3' UTR", 395, c_utr3)
    ],
    'CgiPLT1/2b': [
        ('Promotor', 2000, c_prom), ("5' UTR", 504, c_utr5),
        ('Exón', 158, c_exon), ('Intrón', 174, c_intron),
        ('Exón', 502, c_exon), ('Intrón', 115, c_intron),
        ('Exón', 82, c_exon), ('Intrón', 254, c_intron),
        ('Exón', 11, c_exon), ('Intrón', 128, c_intron),
        ('Exón', 90, c_exon), ('Intrón', 101, c_intron),
        ('Exón', 72, c_exon), ('Intrón', 102, c_intron),
        ('Exón', 130, c_exon), ('Intrón', 117, c_intron),
        ('Exón', 713, c_exon), ("3' UTR", 361, c_utr3)
    ],
    'CgiPLT3/7a': [
        ('Promotor', 2000, c_prom), ("5' UTR", 418, c_utr5),
        ('Exón', 180, c_exon), ('Intrón', 96, c_intron),
        ('Exón', 512, c_exon), ('Intrón', 106, c_intron),
        ('Exón', 84, c_exon), ('Intrón', 1390, c_intron),
        ('Exón', 10, c_exon), ('Intrón', 119, c_intron),
        ('Exón', 90, c_exon), ('Intrón', 119, c_intron),
        ('Exón', 72, c_exon), ('Intrón', 88, c_intron),
        ('Exón', 51, c_exon), ('Intrón', 110, c_intron),
        ('Exón', 78, c_exon), ('Intrón', 211, c_intron),
        ('Exón', 459, c_exon), ("3' UTR", 282, c_utr3)
    ],
    'CgiPLT3/7b': [
        ('Promotor', 2000, c_prom), ("5' UTR", 552, c_utr5),
        ('Exón', 177, c_exon), ('Intrón', 107, c_intron),
        ('Exón', 549, c_exon), ('Intrón', 136, c_intron),
        ('Exón', 85, c_exon), ('Intrón', 1683, c_intron),
        ('Exón', 7, c_exon), ('Intrón', 273, c_intron),
        ('Exón', 91, c_exon), ('Intrón', 123, c_intron),
        ('Exón', 72, c_exon), ('Intrón', 89, c_intron),
        ('Exón', 51, c_exon), ('Intrón', 155, c_intron),
        ('Exón', 78, c_exon), ('Intrón', 108, c_intron),
        ('Exón', 441, c_exon), ("3' UTR", 251, c_utr3)
    ],
    'CgiPLT4': [
        ('Promotor', 2000, c_prom), ("5' UTR", 276, c_utr5),
        ('Exón', 207, c_exon), ('Intrón', 253, c_intron),
        ('Exón', 630, c_exon), ('Intrón', 182, c_intron),
        ('Exón', 84, c_exon), ('Intrón', 88, c_intron),
        ('Exón', 9, c_exon), ('Intrón', 325, c_intron),
        ('Exón', 89, c_exon), ('Intrón', 745, c_intron),
        ('Exón', 73, c_exon), ('Intrón', 95, c_intron),
        ('Exón', 51, c_exon), ('Intrón', 1475, c_intron),
        ('Exón', 79, c_exon), ('Intrón', 96, c_intron),
        ('Exón', 842, c_exon), ("3' UTR", 302, c_utr3)
    ]
}

# ==========================================
# 2. CONFIGURACIÓN Y LECTURA CON ORDEN ESTRICTO
# ==========================================
gen_objetivo = 'CgiPLT3/7b'  # <-- CAMBIA ESTO PARA CADA GEN
archivo_alineamiento = os.path.expanduser(f"~/Downloads/MAFFT_PLT37b_cactus.aln-fasta")

orden_deseado_especies = ['Cgi', 'Ppr', 'Lsc', 'Sth', 'Hun']

try:
    records_raw = list(SeqIO.parse(archivo_alineamiento, "fasta"))
    
    alineamiento = []
    for especie in orden_deseado_especies:
        record = next((r for r in records_raw if especie in r.id), None)
        if record:
            alineamiento.append(record)
        else:
            print(f"Advertencia: No se encontró la especie {especie} en el archivo FASTA.")

    if not alineamiento or 'Cgi' not in alineamiento[0].id:
        raise ValueError("Error crítico: 'Cgi' no se encontró o no quedó en primera posición.")

    molde_record = alineamiento[0]
    seq_molde_coord = str(molde_record.seq) 
    len_ref_real = len(seq_molde_coord.replace('-', ''))
    
    pares = list(itertools.combinations(alineamiento, 2))
    
    print(f"Cargado exitoso. Molde de referencia: {molde_record.id}")
    print(f"Generando {len(pares)} combinaciones para {gen_objetivo}.")
except FileNotFoundError:
    print(f"Error: No se encontró el archivo:\n{archivo_alineamiento}")
    raise

# ==========================================
# 3. CONVERSIÓN DE COORDENADAS (Basadas en Cgi)
# ==========================================
estructura_gen = estructuras_cgi[gen_objetivo]
longitud_upstream = estructura_gen[0][1] + estructura_gen[1][1]
pos_actual = -longitud_upstream
coord_max = -longitud_upstream + len_ref_real
regiones = []

for nombre, longitud, color in estructura_gen:
    inicio = pos_actual
    fin = pos_actual + longitud
    if inicio < coord_max:
        fin_ajustado = min(fin, coord_max)
        regiones.append((inicio, fin_ajustado, nombre, color))
    pos_actual = fin

# ==========================================
# 4. ALGORITMO MULTI-PAR
# ==========================================
def calcular_conservacion_pares_independiente(s_molde, s_A, s_B, inicio_bio, window_size=50):
    """
    Calcula la identidad entre s_A y s_B de forma independiente, 
    usando s_molde (Cgi) únicamente como referencia para las coordenadas del eje X.
    """
    # Usamos las posiciones donde el molde (Cgi) define el sistema de coordenadas,
    # pero evaluamos estrictamente s_A vs s_B en esas columnas del MSA.
    ref_indices = [i for i, base in enumerate(s_molde) if base != '-']
    N_ref = len(ref_indices)
    
    X, Y = [], []
    half_w = window_size // 2
    
    for m in range(N_ref):
        pos_bio = inicio_bio + m
        
        idx_start = max(0, m - half_w)
        idx_end = min(N_ref - 1, m + half_w)
        
        cols_ventana = ref_indices[idx_start : idx_end + 1]

        matches = 0
        total = 0
        for c in cols_ventana:
            a = s_A[c].upper()
            b = s_B[c].upper()
            
            # Si ambas tienen gap en esta columna del MSA global, no aporta a la comparación A-B
            if a == '-' and b == '-':
                continue
                
            total += 1
            # Identidad directa entre A y B
            if a == b and a != 'N':
                matches += 1

        identidad = (matches / total * 100) if total > 0 else 0.0
        X.append(pos_bio)
        Y.append(identidad)
        
    return np.array(X), np.array(Y)

def obtener_etiqueta_par(rec_a, rec_b):
    id_a = rec_a.id
    id_b = rec_b.id
    
    especies = ['Cgi', 'Ppr', 'Lsc', 'Sth', 'Hun']
    code_a = next((sp for sp in especies if sp in id_a), id_a[:3])
    code_b = next((sp for sp in especies if sp in id_b), id_b[:3])
    return f"{code_a}-{code_b}"

# ==========================================
# 5. GENERACIÓN DE LA GRÁFICA (10 Tracks)
# ==========================================
num_tracks = len(pares)

fig, axes = plt.subplots(
    nrows=num_tracks, 
    ncols=1, 
    figsize=(16, 1.4 * num_tracks + 0.6),  # Altura optimizada y más compacta
    sharex=True,
    gridspec_kw={'hspace': 0.05}
)

if num_tracks == 1:
    axes = [axes]

for idx, (ax, (rec_A, rec_B)) in enumerate(zip(axes, pares)):
    seq_A = str(rec_A.seq)
    seq_B = str(rec_B.seq)
    
    X, Y = calcular_conservacion_pares_independiente(seq_molde_coord, seq_A, seq_B, inicio_bio=-longitud_upstream, window_size=50)

    # Relleno de regiones estructurales del molde Cgi
    for inicio, fin, nombre, color in regiones:
        mask = (X >= inicio) & (X <= fin)
        if np.any(mask):
            etiqueta = nombre if idx == 0 else None
            ax.fill_between(X[mask], Y[mask], 0, color=color, alpha=0.85, label=etiqueta)

    # Curva de conservación
    ax.plot(X, Y, color='black', linewidth=0.9, zorder=4)
    
    # Línea del ATG
    etiqueta_atg = 'Inicio AUG (0)' if idx == 0 else None
    ax.axvline(x=0, color='red', linestyle='--', linewidth=1.2, zorder=5, label=etiqueta_atg)

    # Configuración de ejes Y y marcas
    ax.set_ylim(0, 105)
    ax.set_xlim(X[0], X[-1])
    ax.set_yticks([0, 50, 100])
    ax.set_yticklabels(['0%', '50%', '100%'], fontsize=12)

    # Marcas principales del eje X cada 500 pb
    ax.xaxis.set_major_locator(MultipleLocator(500))

    # Etiqueta dinámica e infalible en el lado izquierdo
    etiqueta_izq = obtener_etiqueta_par(rec_A, rec_B)
    ax.set_ylabel(
        etiqueta_izq, 
        fontstyle='italic', 
        fontsize=15, 
        weight='bold', 
        rotation=0, 
        ha='right', 
        va='center', 
        labelpad=15
    )

    ax.grid(True, linestyle=':', alpha=0.3)

# Configuración del eje X en el panel inferior
axes[-1].set_xlabel('Posición relativa al inicio del AUG (pb)', fontsize=15, weight='bold', labelpad=8)

axes[-1].tick_params(axis='x', labelsize=12)  # Cambia '12' al tamaño que prefieras

# Ajustar márgenes para acercar los elementos verticalmente
plt.subplots_adjust(top=0.86, bottom=0.08, left=0.08, right=0.95)

# Título general superior más cerca del borde superior
plt.suptitle(
    f'Conservación a lo largo de las secuencias de los ortólogos {gen_objetivo.replace("Cgi", "")} de los cinco cactus analizados (Ventana = 50 pb)', 
    fontsize=16, 
    weight='bold', 
    y=0.915
)

# ==========================================
# 6. LEYENDA HORIZONTAL SUPERIOR CENTRADA
# ==========================================
handles, labels = axes[0].get_legend_handles_labels()
by_label = dict(zip(labels, handles))

orden_etiquetas = ['Promotor', "5' UTR", 'Exón', 'Intrón', "3' UTR", 'Inicio AUG (0)']
ordered_handles = [by_label[lbl] for lbl in orden_etiquetas if lbl in by_label]
ordered_labels = [lbl for lbl in orden_etiquetas if lbl in by_label]

fig.legend(
    ordered_handles, 
    ordered_labels,
    loc='upper center',
    bbox_to_anchor=(0.50, 0.895),  # Colocado de forma compacta entre el título y las gráficas
    ncol=len(ordered_labels),      
    frameon=False,                # Elimina el rectángulo/caja de la leyenda
    fontsize=16                   # Etiquetas con mayor tamaño y legibilidad
)

# Guardar la imagen
file_path = os.path.join(output_dir, f"mVISTA-like_{gen_objetivo.replace('/', '')}_all_cacti_win50(VER5).jpg")
plt.savefig(file_path, dpi=300, bbox_inches="tight")

plt.show()