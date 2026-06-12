# SimpleRepeatSeq

**Clasificación de secuencias de ADN simples repetidas mediante codificación de k-meros y aprendizaje automático**

Proyecto Integrador · Maestría en Ingeniería Computacional · Bioinformática Computacional

---

## Descripción

SimpleRepeatSeq es un clasificador que, dada una secuencia de ADN, determina si es un *simple repeat* (telómeros, expansiones de enfermedades como Huntington/X-Frágil/ELA, STRs forenses, satélites centroméricos) o una secuencia genómica normal (exones, intrones, promotores, UTRs, etc.).

El sistema evalúa sistemáticamente **6 encodings × 3 preprocesadores × 8 modelos** mediante validación cruzada estratificada 10-fold.

---

## Estructura del repositorio

```
SimpleRepeatSeq/
├── data/
│   ├── SimpleRepeatSeq_dataset.fasta    ← coloca el dataset aquí
│   └── encoded/                         ← arrays .npy generados por notebook 01
├── notebooks/
│   ├── 01_load_data.ipynb              ← Carga y generación de encodings
│   └── 02_experiments.ipynb            ← Experimentos, métricas y figuras
├── figures/                             ← Figuras generadas automáticamente
├── models/
│   └── best_model.pkl                   ← Mejor modelo entrenado (RF + kmer4 + scaled)
├── article/
│   └── paper.pdf                        ← Artículo científico
├── classify.py                          ← CLI de clasificación
├── requirements.txt
└── README.md
```

---

## Instalación local (entorno virtual)

**Requisitos previos:** Python ≥ 3.9, `git`

> **Debian/Ubuntu:** si `python3 -m venv` falla con *"ensurepip is not available"*, instala el módulo venv antes de continuar:
> ```bash
> sudo apt install python3.12-venv -y
> ```

```bash
# 1. Clonar el repositorio
git clone https://github.com/leandropa00/SimpleRepeatSeq.git
cd SimpleRepeatSeq

# 2. Crear y activar el entorno virtual
python3 -m venv .venv
source .venv/bin/activate          # Linux / macOS
# .venv\Scripts\activate           # Windows (PowerShell)

# 3. Instalar dependencias
pip install --upgrade pip
pip install -r requirements.txt

# 4. Registrar el kernel de Jupyter para este entorno
python -m ipykernel install --user --name simplerepeatseq --display-name "SimpleRepeatSeq"
```

**Dependencias:** biopython, numpy, pandas, matplotlib, seaborn, scikit-learn, xgboost, tensorflow ≥ 2.13, scikeras, joblib.

---

## Uso — Notebooks

```bash
# Con el entorno virtual activo:
jupyter notebook
```

Abrir en el navegador y ejecutar en orden:

1. `notebooks/01_load_data.ipynb` — carga el FASTA y genera los 6 encodings (arrays `.npy` en `data/encoded/`)
2. `notebooks/02_experiments.ipynb` — experimentos completos (6 encodings × 3 preprocesadores × 8 modelos), figuras y modelo final

> **Nota:** al abrir cada notebook, seleccionar el kernel **SimpleRepeatSeq** (menú Kernel → Change kernel).  
> El notebook 02 tarda aproximadamente **15–40 min** dependiendo del hardware.

Los notebooks también son compatibles con **Google Colab**: clonar el repo, subir el FASTA a `data/` y ejecutar las celdas de instalación incluidas al inicio.

---

## Uso — CLI de clasificación

```bash
python classify.py --fasta data/SimpleRepeatSeq_dataset.fasta \
                   --model models/best_model.pkl \
                   --encoding kmer4
```

**Parámetros:**

| Parámetro | Descripción |
|-----------|-------------|
| `--fasta` | Archivo FASTA de entrada |
| `--model` | Ruta al modelo `.pkl` entrenado |
| `--encoding` | Esquema de codificación: `kmer3`, `kmer4`, `kmer6`, `dax`, `eiip`, `complementary` (o alias `km3`, `km4`, `km6`) |
| `--k` | Alternativa a `--encoding` para k-meros; valores válidos: 3, 4, 6 |

**Salida de ejemplo:**

```
ID                                                 Clase predicha       Probabilidad
-------------------------------------------------------------------------------------
CAG_huntington_007#simpleSeq#210                   simpleSeq                  1.0000
exon_coding_042#no_simpleSeq#175                   no_simpleSeq               1.0000
```

---

## Resultados principales

Validación cruzada estratificada 10-fold (888 secuencias, dataset balanceado).

| Encoding | Mejor modelo | Preproc | Accuracy | AUC   |
|----------|-------------|---------|----------|-------|
| kmer4    | RF          | scaled  | 1.000    | 1.000 |
| kmer3    | RF          | raw     | 1.000    | 1.000 |
| kmer6    | SVM         | pca     | 1.000    | 1.000 |
| eiip     | NB          | pca     | 0.985    | 1.000 |
| dax      | XGB         | raw     | 0.987    | 0.999 |
| complementary | SVM  | raw     | 0.989    | 1.000 |

**Mejor configuración global:** Random Forest + k-meros k=4 + StandardScaler  
Accuracy = 1.000 · AUC = 1.000 · F1-weighted = 1.000

**Ranking de encodings por AUC medio:** kmer4 > kmer3 > kmer6 > eiip > dax > complementary  
**Ranking de modelos por AUC medio:** RF > XGB > SVM > LR > LDA > NB > KNN > DNN

---

## Dataset

888 secuencias balanceadas (444 *simpleSeq* + 444 *no_simpleSeq*), longitudes 100–500 pb con 4–12% de ruido de secuenciación.

- **simpleSeq:** telómeros (TTAGGG), expansiones de enfermedades (CAG/Huntington, CGG/X-Frágil, GGGGCC/ELA), STRs forenses, satélites centroméricos, microsatélites di/trinucleotídicos.
- **no_simpleSeq:** exones, intrones, promotores CpG, 5'UTR, 3'UTR, enhancers, mitocondrial, SINEs.

Formato de identificador: `>CAG_huntington_007#simpleSeq#210`
