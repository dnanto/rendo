#!/usr/bin/env python3

import re
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
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


def digest(enzyme, c5, c3, string):
	positions = [ele.start() for ele in re.finditer(enzyme, string)]

	s5 = [0] + [pos + c5 for pos in positions] + [len(string) + 1]
	s3 = [0] + [pos + c3 for pos in positions] + [len(string) + 1]

	for idx in range(len(s5) - 1):
		yield string[s5[idx]:s5[idx + 1]], Seq.complement(string[s3[idx]:s3[idx + 1]])


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

	enzyme = args.enzyme.upper()
	decoded = "".join(decode(enzyme.replace("<>", "|").replace("><", "|")))
	regex = "".join(regexify(enzyme))

	print("0123456789" * int(len(decoded) / 10 + 1))
	print(decoded)
	print(regex)

	assert re.match(regex, decoded.replace(">", "").replace("<", "").replace("|", ""))

	n1, n2, n3 = (decoded.count(ch) for ch in "|><")
	if n1:
		if n2 or n3:
			raise ValueError("enzyme definition error...")
	elif not (n2 == 1 and n3 == 1):
		raise ValueError("enzyme definition error...")

	print("blunt" if n1 else "sticky")

	c5, c3 = "||" if n1 else "><"
	c5, c3 = decoded.index(c5), decoded.index(c3)
	c5 -= c3 < c5
	c3 -= c5 < c3

	print("5' cut @", c5)
	print("3' cut @", c3)

	with args.file as file:
		for record in SeqIO.parse(file, args.fmt):
			print(record)
			r5 = []
			r3 = []
			for s5, s3 in digest(regex, c5, c3, str(record.seq)):
				print(s5, s3)
				r5.append(s5)
				r3.append(s3)
			print("".join(r5))
			print("".join(r3))
			assert "".join(r5) == str(record.seq)
			assert "".join(r3) == Seq.complement(str(record.seq))

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
