#!/usr/bin/env python3

import re
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from signal import signal, SIGPIPE, SIG_DFL

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


def parse_argv(argv):
	parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

	parser.add_argument(
		"file", type=FileType()
	)
	parser.add_argument(
		"enzyme"
	)
	parser.add_argument(
		"-fmt", "--fmt", default="fasta"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	# with args.file as file:
	# 	for record in SeqIO.parse(file, args.fmt):
	# 		print(record.id)

	enzyme = args.enzyme.upper()
	decoded = "".join(decode(enzyme.replace("<>", "|").replace("><", "|")))
	regex = "".join(regexify(enzyme))

	print("0123456789" * int(len(decoded) / 10 + 1))
	print(decoded)
	print(regex)

	assert re.match(regex, decoded.replace(">", "").replace("<", "").replace("|", ""))

	n1, n2, n3 = (decoded.count(ch) for ch in "|><")
	if (n1 and n1 != 1 or n2 or n3) or not n1 and n2 != 1 and n3 != 1:
		raise ValueError("enzyme definition error...")

	print("blunt" if n1 else "sticky")

	c5, c3 = "||" if n1 else "><"
	c5, c3 = decoded.index(c5), decoded.index(c3)

	print("5' cut @ , ", c5)
	print("3' cut @ , ", c3)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
