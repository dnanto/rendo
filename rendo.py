#!/usr/bin/env python3

import re
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from signal import signal, SIGPIPE, SIG_DFL

from Bio import Restriction
from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import Seq


def decode(enzyme):
	for ele in re.finditer(r"\w\d+|\w|[^\w\d]+", enzyme):
		token = ele.group(0)
		base, rep = token[0], token[1:]
		yield base * (int(rep) if rep else 1)


def regexify(enzyme):
	for ele in re.finditer(r"\w\d+|\w", enzyme):
		token = ele.group(0)
		base, rep = token[0], token[1:]

		bases = ambiguous_dna_values[base]
		bases = f"[{bases + base}]" if len(bases) > 1 else bases

		rep = f"{{{int(rep)}}}" if rep else rep

		yield bases + rep


def load_enzyme(enzyme):
	enzyme = enzyme.upper()
	decoded = "".join(decode(enzyme.replace("<>", "|").replace("><", "|")))
	regex = "".join(regexify(enzyme))

	assert re.match(regex, decoded.replace(">", "").replace("<", "").replace("|", ""))

	n1, n2, n3 = (decoded.count(ch) for ch in "|><")
	if n1:
		if n2 or n3:
			raise ValueError("enzyme definition error...")
	elif not (n2 == 1 and n3 == 1):
		raise ValueError("enzyme definition error...")

	print(decoded, "->", regex, "(blunt)" if n1 else "(sticky)", file=sys.stderr)

	c5, c3 = "||" if n1 else "><"
	c5, c3 = decoded.index(c5), decoded.index(c3)
	c5 -= c3 < c5
	c3 -= c5 < c3

	print("5' cut @", c5, file=sys.stderr)
	print("3' cut @", c3, file=sys.stderr)

	return regex, c5, c3


def digest(enzyme, c5, c3, dna, fragments=True):
	positions = [ele.start() for ele in re.finditer(enzyme, dna)]
	p5 = [0] + [pos + c5 for pos in positions] + [len(dna) + 1]
	p3 = [0] + [pos + c3 for pos in positions] + [len(dna) + 1]
	indexes = range(len(p5) - 1)
	if fragments:
		for idx in indexes:
			s5, s3 = dna[p5[idx]:p5[idx + 1]], Seq.complement(dna[p3[idx]:p3[idx + 1]])
			yield p5[idx], p5[idx + 1], p3[idx], p3[idx + 1], s5, s3
	else:
		for idx in indexes:
			l5, l3 = p5[idx + 1] - p5[idx] + 1, p3[idx + 1] - p3[idx] + 1
			yield p5[idx], p5[idx + 1], p3[idx], p3[idx + 1], l5, l3


def process(argv):
	args = parse_argv(argv[1:])

	with args.enzyme.open() as file:
		enzymes = list((tag, *load_enzyme(enzyme)) for tag, enzyme in map(str.split, file))

	with args.file as file1, args.out as file2:
		print("id", "enzyme", "regex", "c1p5", "c2p5", "c1p3", "c2p3", "r5", "r3", sep="\t", file=file2)
		for record in SeqIO.parse(file1, args.fmt):
			for tag, regex, c5, c3 in enzymes:
				print(record.id, ">", tag, file=sys.stderr)
				results = digest(regex, c5, c3, str(record.seq), args.fragments)
				results = ((record.id, tag, regex, *row) for row in results)
				writer(file2, delimiter="\t").writerows(results)


def parse_argv(argv):
	topology = ("linear", "circular")
	parser = ArgumentParser(description="calculate restriction fragment sizes", formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("file", type=FileType(), help="the sequence file")
	parser.add_argument("-enzymes", nargs="+", default=["EcoRI"], help="the list of enzyme names")
	parser.add_argument("-topology", choices=topology, default=topology[0], help="the sequence topology")
	args = parser.parse_args(argv)
	return args


def main(argv):
	args = parse_argv(argv[1:])

	enzymes = [getattr(Restriction, ele) for ele in args.enzymes]
	is_linear = args.topology == "linear"

	with args.file as file:
		print("id", "enzyme", "size", sep="\t")
		for record in SeqIO.parse(file, "fasta"):
			seq = Seq(str(record.seq).replace("?", "N"), alphabet=IUPACAmbiguousDNA())
			for enzyme in enzymes:
				for frag in enzyme.catalyze(seq, linear=is_linear):
					print(record.id, enzyme, len(frag), sep="\t")

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
