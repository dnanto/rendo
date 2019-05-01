#!/usr/bin/env python3

import re
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from csv import writer
from pathlib import Path
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO, Seq
from Bio.Data.IUPACData import ambiguous_dna_values


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
	p5 = [0] + [pos + c5 for pos in positions] + [len(dna)]
	p3 = [0] + [pos + c3 for pos in positions] + [len(dna)]
	indexes = range(len(p5) - 1)
	if fragments:
		for idx in indexes:
			s5, s3 = dna[p5[idx]:p5[idx + 1]], Seq.complement(dna[p3[idx]:p3[idx + 1]])
			yield p5[idx] + 1, p5[idx + 1], p3[idx], p3[idx + 1], s5, s3
	else:
		for idx in indexes:
			l5, l3 = p5[idx + 1] - p5[idx] + 1, p3[idx + 1] - p3[idx] + 1
			yield p5[idx] + 1, p5[idx + 1], p3[idx], p3[idx + 1], l5, l3


def parse_argv(argv):
	parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

	parser.add_argument(
		"file", type=FileType()
	)
	parser.add_argument(
		"enzyme", type=Path
	)
	parser.add_argument(
		"-fmt", "--fmt", default="fasta"
	)
	parser.add_argument(
		"-from-file", "--from-file", action="store_true"
	)
	parser.add_argument(
		"-out", "--out", type=FileType("w"), default="-"
	)
	parser.add_argument(
		"-fmt-out", "--fmt-out", default="tab"
	)
	parser.add_argument(
		"-fragments", "--fragments", action="store_true"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	with args.enzyme.open() as file:
		enzymes = (line for line in map(str.strip, file) if line[0] != "#")
		enzymes = [(key, *load_enzyme(val)) for key, val in map(str.split, enzymes)]

	with args.file as file1, args.out as file2:
		print("id", "enzyme", "regex", "c1p5", "c2p5", "c1p3", "c2p3", "r5", "r3", sep="\t", file=file2)
		for record in SeqIO.parse(file1, args.fmt):
			for tag, regex, c5, c3 in enzymes:
				print(record.id, ">", tag, file=sys.stderr)
				results = digest(regex, c5, c3, str(record.seq), args.fragments)
				results = ((record.id, tag, regex, *row) for row in results)
				writer(file2, delimiter="\t").writerows(results)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
