cat pangenome_mapping | grep + | cut -d, -f1,4 | sed -E 's/("([^"]*)")?,/\2\t/g' >> parsed_pangenome_mapping


