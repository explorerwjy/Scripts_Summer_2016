#import vcf.filters
import vcf.filters

class ExAC_NFE_MAF(vcf.filters.Base):
    "Filter sites by ExAC NFE MAF value"
    
    name = "exac-nfe-maf"
    
    @classmethod
    def customize_parser(self, parser):
        parser.add_argument("--site-ExAC-NFE-MAF", type = int, default = 0.01, help = "Filter sites above this MAF value")
        
    def __init__(self, args):
        self.threshold = args.site_ExAC_NFE_MAF
        
    def __call__(self, record):
        for maf in record.INFO["ExAC_NFE"]:
            if float(0 if maf is None else maf) > self.threshold:
                return maf
                
class Func_refGene(vcf.filters.Base):
    "Function of gene can be exonic"
    
    name = "func-refgene-exonic"
    
    @classmethod
    def customize_parser(self, parser):
        parser.add_argument("--func-refGene-exonic", type = str, default = "exonic", help = "Filter sites for predicted function of gene to be exonic")
        
    def __init__(self, args):
        self.threshold = args.func_refGene
        
    def __call__(self, record):
        for gene_function in record.INFO["Func.refGene"]:
            if gene_function == "exonic":
                return gene_function

class Func_refGene(vcf.filters.Base):
    "Function of gene can be splicing"
    
    name = "func-refgene-splicing"
    
    @classmethod
    def customize_parser(self, parser):
        parser.add_argument("--func-refGene-splicing", type = str, default = "exonic", help = "Filter sites for predicted function of gene to be exonic")
        
    def __init__(self, args):
        self.threshold = args.func_refGene
        
    def __call__(self, record):
        for gene_function in record.INFO["Func.refGene"]:
            if gene_function == "exonic":
                return gene_function
class Func_refGene(vcf.filters.Base):
    "Filter sites by predicted function of gene"
    
    name = "func-refgene-exonic"
    
    @classmethod
    def customize_parser(self, parser):
        parser.add_argument("--func-refGene", type = str, default = "exonic", help = "Filter sites for predicted function of gene to be exonic")
        
    def __init__(self, args):
        self.threshold = args.func_refGene
        
    def __call__(self, record):
        for gene_function in record.INFO["Func.refGene"]:
            if gene_function == "exonic":
                return gene_function
class ExonicFunc_refGene(vcf.filters.Base):
    "Filter sites by predicted function of gene"
    
    name = "exonicfunc-refgene"
    
    @classmethod
    def customize_parser(self, parser):
        parser.add_argument("--exonicfunc-refGene", type = str, default = "nonsynonymous_SNV", help = "Filter sites for predicted function of gene to be exonic")
        
    def __init__(self, args):
        self.threshold = args.exonicfunc_refGene
        
    def __call__(self, record):
        for exonicgene_function in record.INFO["ExonicFunc.refGene"]:
            if str("exclude" if exonicgene_function is None else exonicgene_function) != self.threshold:
                return str("exclude" if exonicgene_function is None else exonicgene_function)

class Not_synonymous(vcf.filters.Base):
	"Filter sites by predicted function of gene"

	name = "exonicfunc-refGene-not-syn"

	@classmethod
	def customize_parser(self, parser):
		parser.add_argument("--exonicfunc-refGene-not-syn", type = str, default = "synonymous_SNV", help = "Filter sites for predicted function of gene to be exonic")

	def __init__(self, args):
		self.threshold = args.exonicfunc_refGene_not_syn
	
	def __call__(self, record):
		for exonicgene_function in record.INFO["ExonicFunc.refGene"]:
			if str("exclude" if exonicgene_function is None else exonicgene_function) == self.threshold:
				return str("exclude" if exonicgene_function is None else exonicgene_function)
