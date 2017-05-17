LOAD CSV WITH HEADERS FROM "file:///cpg.csv" AS csvLine 
CREATE (t1:cpg { 
	name: csvLine.cpg, 
	chr: csvLine.chr, 
	pos: toInteger(csvLine.pos)
})

LOAD CSV WITH HEADERS FROM "file:///snp.csv" AS csvLine 
CREATE (t1:snp { 
	name: csvLine.snp, 
	chr: csvLine.chr, 
	pos: toInteger(csvLine.pos), 
	maf: toFloat(csvLine.maf), 
	type: csvLine.type 
})

LOAD CSV WITH HEADERS FROM "file:///snp_cpg.csv" AS csvLine 
MATCH 
(t1:snp { 
	name: csvLine.snp }),
(t2:cpg { 
	name: csvLine.cpg }) 
CREATE (t1)-[m:GA { 
	b: toFloat(csvLine.b), 
	pval: toFloat(csvLine.pval) 
}]->(t2)

LOAD CSV WITH HEADERS FROM "file:///cpg_cpg.csv" AS csvLine 
MATCH 
(t1:cpg { 
	name: csvLine.cpg1 }),
(t2:cpg { 
	name: csvLine.cpg2 }) 
CREATE (t1)-[m:MR { 
	method: csvLine.method,
	b: toFloat(csvLine.waldratio), 
	sign: toInteger(csvLine.pval)
}]->(t2)

CREATE INDEX ON :cpg
CREATE INDEX ON :snp



~/Downloads/neo4j-community-3.2.0/bin/neo4j-admin import \
	--database cis-trans.db \
	--id-type string \
	--nodes:CpG ~/repo/godmc_phase2_analysis/05_cis-trans-networks/data/cpg.csv \
	--nodes:SNP ~/repo/godmc_phase2_analysis/05_cis-trans-networks/data/snp.csv \
	--relationships:GA ~/repo/godmc_phase2_analysis/05_cis-trans-networks/data/snp_cpg.csv \
	--relationships:MR ~/repo/godmc_phase2_analysis/05_cis-trans-networks/data/cpg_cpg.csv
