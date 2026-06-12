#!/usr/bin/env python3
"""SimpleRepeatSeq — CLI de clasificación.

Clasifica cada secuencia de un archivo FASTA usando el mejor modelo entrenado.
Imprime: ID de secuencia, clase predicha y probabilidad.

Uso:
    python classify.py --fasta input.fasta --model models/best_model.pkl --encoding kmer4
    python classify.py --fasta input.fasta --model models/best_model.pkl --k 4
"""
import argparse
import sys
from itertools import product

import joblib
import numpy as np
from Bio import SeqIO


ENCODING_ALIASES = {
    'km3': 'kmer3', 'km4': 'kmer4', 'km6': 'kmer6',
    'kmer3': 'kmer3', 'kmer4': 'kmer4', 'kmer6': 'kmer6',
    'dax': 'dax', 'eiip': 'eiip', 'complementary': 'complementary',
}


def _pad(vec, max_len, pad_value=0.0):
    if len(vec) >= max_len:
        return vec[:max_len]
    return vec + [pad_value] * (max_len - len(vec))


def _normalizar(sequence):
    return ''.join(nt if nt in 'ACGT' else 'N' for nt in sequence.upper())


def dax_encode(sequence, max_len):
    mapping = {'A': 2.0, 'C': 0.0, 'G': 3.0, 'T': 1.0}
    return _pad([mapping.get(nt, 1.5) for nt in sequence.upper()], max_len)


def eiip_encode(sequence, max_len):
    mapping = {'A': 0.126, 'C': 0.134, 'G': 0.081, 'T': 0.134}
    return _pad([mapping.get(nt, 0.0) for nt in sequence.upper()], max_len)


def complementary_encode(sequence, max_len):
    mapping = {'A': 2.0, 'C': -1.0, 'G': 1.0, 'T': -2.0}
    return _pad([mapping.get(nt, 0.0) for nt in sequence.upper()], max_len)


def kmer_encode(sequence, k):
    alfabeto = 'ACGTN'
    vocab = sorted(''.join(p) for p in product(alfabeto, repeat=k))
    kmer_to_idx = {km: i for i, km in enumerate(vocab)}
    counts = np.zeros(len(vocab), dtype=np.float32)
    seq = _normalizar(sequence)
    n_windows = len(seq) - k + 1
    if n_windows > 0:
        for j in range(n_windows):
            km = seq[j:j+k]
            if km in kmer_to_idx:
                counts[kmer_to_idx[km]] += 1
        counts /= n_windows
    return counts.tolist()


def encode_sequence(sequence, encoding, max_len):
    seq = str(sequence)
    if encoding == 'dax':
        return dax_encode(seq, max_len)
    if encoding == 'eiip':
        return eiip_encode(seq, max_len)
    if encoding == 'complementary':
        return complementary_encode(seq, max_len)
    if encoding == 'kmer3':
        return kmer_encode(seq, 3)
    if encoding == 'kmer4':
        return kmer_encode(seq, 4)
    if encoding == 'kmer6':
        return kmer_encode(seq, 6)
    raise ValueError(f'Encoding desconocido: {encoding}')


def main():
    parser = argparse.ArgumentParser(
        description='SimpleRepeatSeq: clasifica secuencias de ADN de un archivo FASTA.'
    )
    parser.add_argument('--fasta',     required=True,
                        help='Archivo FASTA de entrada')
    parser.add_argument('--model',     required=True,
                        help='Ruta al modelo entrenado (.pkl)')
    parser.add_argument('--encoding',  default=None,
                        help='Esquema de codificación: kmer3 | kmer4 | kmer6 | '
                             'dax | eiip | complementary  (o alias: km3 | km4 | km6)')
    parser.add_argument('--k',         type=int, default=None,
                        help='Valor de k para k-meros (alternativa a --encoding; '
                             'valores válidos: 3, 4, 6)')
    args = parser.parse_args()

    # Resolver encoding
    if args.encoding is not None:
        enc = ENCODING_ALIASES.get(args.encoding.lower())
        if enc is None:
            print(f'ERROR: encoding desconocido: "{args.encoding}".\n'
                  f'Valores válidos: {", ".join(ENCODING_ALIASES.keys())}',
                  file=sys.stderr)
            sys.exit(1)
    elif args.k is not None:
        enc = f'kmer{args.k}'
        if enc not in ('kmer3', 'kmer4', 'kmer6'):
            print(f'ERROR: --k {args.k} no soportado. Valores válidos: 3, 4, 6',
                  file=sys.stderr)
            sys.exit(1)
    else:
        print('ERROR: especifica --encoding o --k.', file=sys.stderr)
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Cargar modelo
    try:
        bundle   = joblib.load(args.model)
        pipeline = bundle['pipeline']
        meta     = bundle['meta']
    except Exception as e:
        print(f'ERROR al cargar el modelo: {e}', file=sys.stderr)
        sys.exit(1)

    max_len = meta.get('max_len', 500)
    clases  = meta.get('clases', ['no_simpleSeq', 'simpleSeq'])
    model_enc = meta.get('encoding', '')

    if model_enc and model_enc != enc:
        print(f'AVISO: el modelo fue entrenado con encoding="{model_enc}" '
              f'pero se usa encoding="{enc}". Los resultados pueden ser incorrectos.',
              file=sys.stderr)

    # Cargar secuencias
    try:
        records = list(SeqIO.parse(args.fasta, 'fasta'))
    except Exception as e:
        print(f'ERROR al leer el FASTA: {e}', file=sys.stderr)
        sys.exit(1)

    if not records:
        print('ERROR: no se encontraron secuencias en el archivo FASTA.', file=sys.stderr)
        sys.exit(1)

    # Clasificar y mostrar resultados
    header = f'{"ID":<50} {"Clase predicha":<20} {"Probabilidad":>12}'
    print(header)
    print('-' * len(header))

    for rec in records:
        feat = encode_sequence(rec.seq, enc, max_len)
        X    = np.array([feat], dtype=np.float32)
        pred      = int(pipeline.predict(X)[0])
        proba_all = pipeline.predict_proba(X)[0]
        clase     = clases[pred]
        prob_val  = proba_all[pred]
        print(f'{rec.id:<50} {clase:<20} {prob_val:>12.4f}')


if __name__ == '__main__':
    main()
