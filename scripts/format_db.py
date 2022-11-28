# -------------------------------
# Title: format_db.py
# Author: Silver A. Wolf
# Last Modified: Mon, 28.11.2022
# Version: 0.0.1
# -------------------------------

new_db = open("db/megares_filtered.fasta", "w")
with open("db/megares_full_database_v2.00.fasta") as db:
	for entry in db:
		if entry[0] == ">":
			if "RequiresSNPConfirmation" not in entry:
				reading = True
				new_db.write(entry)
			else:
				reading = False
		else:
			if reading == True:
				new_db.write(entry)
new_db.close()

new_csv = open("db/megares_filtered.csv", "w")
with open("db/megares_drugs_annotations_v2.00.csv") as csv:
	for line in csv:
		if "RequiresSNPConfirmation" not in line:
			new_csv.write(line)
new_csv.close()