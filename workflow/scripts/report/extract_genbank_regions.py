#!/usr/bin/env python3

import logging
import json

from Bio import SeqIO
from Bio.SeqFeature import ExactPosition


if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"],
        level=logging.INFO
    )

    logging.info("Reading GenBank file")
    gb = SeqIO.read(snakemake.input.gb, format="gb")

    logging.info("Calculating CDS regions")
    regions = {}
    for feature in gb.features:
        identifier = "|".join(feature.qualifiers.get(snakemake.params.gb_qualifier, []))
        if identifier == "":
            logging.error(f"Feature at {feature.location} has no qualifier '{snakemake.params.gb_qualifier_display}'")
        elif identifier in regions:
            logging.warning(f"Identifier '{identifier}' is already among coding records and will not be replaced by feature at {feature.location}")
        else:
            logging.debug(f"Adding feature")
            if type(feature.location.start) is not ExactPosition or type(feature.location.start) is not ExactPosition:
                logging.warning(f"Feature '{identifier}' location is not exact but will be treated as such: {feature.location}")
            regions[identifier] = (
                int(feature.location.start.real + 1),
                int(feature.location.end.real)
            )
    
    logging.info("Calculating intergenic regions")
    cds_names = tuple(regions.keys())
    intergenic_count = 0
    for i in range(1, len(cds_names)):
        start1, end1 = regions[cds_names[i-1]]
        start2, end2 = regions[cds_names[i]]
        if (start2 - end1) > 1:
            intergenic_count += 1
            regions[f"Intergenic_{intergenic_count}"] = (end1 + 1, start2 - 1)
    
    logging.info("Writing regions to JSON file")
    with open(snakemake.output.regions, "w") as fw:
        json.dump(regions, fw, indent=2)
