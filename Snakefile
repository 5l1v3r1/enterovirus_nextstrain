
rule files:
    params:
        dropped_strains = "config/dropped_strains.txt",
        reference = "config/ev_d68_reference_genome.gb",
        colors = "config/colors.tsv",
        auspice_config = "config/auspice_config.json"

files = rules.files.params

rule concat_meta:
    input:
        metadata = ["data/EVD68_from_sweden.csv", "data/EVD68_minor_from_sweden.csv", "data/EVD68_VIPR.csv"]
    output:
        metadata = "results/metadata.tsv"
    run:
        import pandas as pd
        from augur.parse import fix_dates, forbidden_characters
        md = []
        for fname in input.metadata:
            tmp = pd.read_csv(fname, sep='\t' if fname.endswith('tsv') else ',')
            tmp_dates = []
            for x in tmp.date:
                tmp_dates.append(fix_dates(x, dayfirst='VIPR' not  in fname))
            tmp.date = tmp_dates
            tmp_name = []
            for x in tmp.strain:
                f = x
                for c,r in forbidden_characters:
                    f=f.replace(c,r)
                tmp_name.append(f)
            tmp.strain = tmp_name
            md.append(tmp)
        all_meta = pd.concat(md)
        all_meta.to_csv(output.metadata, sep='\t')

rule concat_seqs:
    input:
        fasta = ["data/ev_d68_genomes.fasta", "data/ev_d68_genomes_sweden.fasta", "data/ev_d68_minor_genomes_sweden.fasta"]
    output:
        sequences = "results/sequences.fasta"
    shell:
        """
        cat {input.fasta} > {output.sequences}
        """

rule filter:
    input:
        sequences = rules.concat_seqs.output.sequences,
        metadata = rules.concat_meta.output.metadata,
        exclude = files.dropped_strains
    output:
        sequences = "results/filtered.fasta"
    params:
        sequences_per_category = 200,
        categories = "country year month",
        min_date = 2000
    shell:
        """
        augur filter --sequences {input.sequences} --metadata {input.metadata} \
            --output {output.sequences} \
            --group-by {params.categories} \
            --sequences-per-group {params.sequences_per_category} \
            --exclude {input.exclude}  --min-date {params.min_date}
        """

rule align:
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned_genome.fasta"
    shell:
        """
        augur align --sequences {input.sequences} --output {output.alignment} \
            --reference-sequence {input.reference} --fill-gaps
        """

rule sub_alignments:
    input:
        rules.align.output
    output:
        subaln = ["results/aligned_vp4vp1.fasta", "results/aligned_p2p3.fasta"]
    params:
        boundaries = [(699,3282), (3282,7263)]
    run:
        from Bio import AlignIO
        aln = AlignIO.read(input[0], 'fasta')
        print(aln)
        for fname,b in zip(output.subaln, params.boundaries):
            print(b)
            AlignIO.write(aln[:,b[0]:b[1]], fname, 'fasta')


rule tree:
    input:
        alignment = [rules.align.output.alignment] + rules.sub_alignments.output
    output:
        tree = "results/raw_tree_{seg}.nwk"
    shell:
        """
        augur tree --alignment results/aligned_{wildcards.seg}.fasta --output {output.tree}
        """

rule refine:
    input:
        tree = rules.tree.output.tree,
        alignment = 'results/aligned_{seg}.fasta',
        metadata = rules.concat_meta.output.metadata
    output:
        tree = "results/tree_{seg}.nwk",
        node_data = "results/branch_lengths_{seg}.json"
    params:
        clock_filter_iqd = 4
    shell:
        """
        augur refine --tree {input.tree} --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} --output-node-data {output.node_data} \
            --timetree --date-confidence --date-inference marginal --coalescent opt \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    input:
        tree = rules.refine.output.tree,
        alignment = 'results/aligned_{seg}.fasta',
    output:
        node_data = "results/nt_muts_{seg}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral --tree {input.tree} --alignment {input.alignment} \
            --output {output.node_data} --inference {params.inference}
        """

rule translate:
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = 'config/ev_d68_reference_{seg}.gb'
    output:
        node_data = "results/aa_muts_{seg}.json"
    shell:
        """
        augur translate --tree {input.tree} --ancestral-sequences {input.node_data} \
            --output {output.node_data} --reference-sequence {input.reference}
        """

rule traits:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.concat_meta.output.metadata
    output:
        node_data = "results/traits_{seg}.json",
    params:
        columns = "country"
    shell:
        """
        augur traits --tree {input.tree} --metadata {input.metadata} \
            --output {output.node_data} --confidence --columns {params.columns}
        """

rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.concat_meta.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = files.colors,
        auspice_config = files.auspice_config
    output:
        auspice_tree = "auspice/ev_d68_{seg}_tree.json",
        auspice_meta = "auspice/ev_d68_{seg}_meta.json"
    shell:
        """
        augur export --tree {input.tree} --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --colors {input.colors} --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} --output-meta {output.auspice_meta}
        """
