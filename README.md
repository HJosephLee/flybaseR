# to install

```
devtools::install_github('hangnoh/flybaseR')
```


# id.converter
Converting or Updating FlyBase IDs or Symbols

Description

The function takes FlyBase IDs (e.g. FBgn0000003) or Gene symbols as an input, and converts it into updated IDs or symbols using the FlyBase ID converter (web). The function accesses FlyBase, so requires internet-connection. FlyBase ID inputs are bundled as 1,000. 100 for symbols. FlyBase IDs for genes that are split into multiple genes will be concatenated with two colons (::). Genes that does not have matching IDs will be shown as "unknown". Certain gene symbols would appear as "unknown" even if the gene exists, and have FlyBase IDs. This is because the ID converter in FlyBase website cannot convert the gene. For example, CG31976 cannot be converted by FlyBase, although you can find the gene from the gene report. Setting diehard.symbols = T will look for gene report pages of such unconvertible genes one by one. The process is essentially slow because it accesses FlyBase for each gene.

Usage
```
id.converter(x, symbols = F, bundle.size = 1000, DmelOnly = T, polite.access = 0, diehard.symbols = F, convert.into)
Arguments
```
x
> a vector. FlyBase IDs or names to be converted.
symbols
> Logical. If TRUE, the output will be gene symbols, rather than FlyBase IDs. Default = F
bundle.size
> Numeric. The number of FlyBase IDs or symbols to be submitted to FlyBase at once. Default is 1,000 if there are less than 100 symbols; 100 if more than 1,000 symbols. Reduce the number down if Timeout error occurs.
DmelOnly
> Logical. If TRUE, non-melanogaster gene IDs will be ignored. Default = T.
polite.access
> Numeric. Intervals between FlyBase access for each bundle as seconds. Default = 0.
diehard.symbols
> Logical. If TRUE, ntervals between FlyBase access for each bundle as seconds. Default = 0.
convert.into
> "genes", "transcripts", or "polypeptides". "g", "t", or "p" is also possible. If missing, the IDs will be updated to the most recent IDs only.

Examples
```
id.converter(x, symbols = T)
id.converter(x, bundle.size = 50, be.polite = 10, convert.into = "transcripts")
id.converter(x, symbols = T, bundle.size = 50, diehard.symbols = T)
```


# id.converter2
Updating FlyBase IDs to a certain version.

Description

The function takes FlyBase IDs (e.g. FBgn#######) as an input, and converts it into certain versions of IDs. This function is not able to handle gene symbols. The function accesses the FlyBase FTP site, so requires internet-connection. FlyBase IDs for genes that are split into multiple genes will be concatenated with two colons (::). Genes that does not have matching IDs will be shown as "unknown".

Usage
```
id.converter2(x, version, thread = 1)
Arguments

x- a vector. FlyBase IDs.
version - FlyBase ID version for the updated result. It should be either version numbers (e.g. 6.12) or FlyBase Release Numbers (e.g. FB2017_04). The default is "current".
thread - The function itself is slow, so in order to speed up you can use multiple CPU threads for parallelization. Default : thread = 1.
```

Examples
```
id.converter2(x, version=6.12)
id.converter2(x, version="FB2016_04", thread=4)
```


# j2g
Open FlyBase gene reports using a web browser

Description

The function takes FlyBase IDs (e.g. FBgn#######) or Gene symbols as an input, and opens FlyBase gene report from the default web browser (e.g. Jump2Gene).

Usage
```
j2g(x, n = 10)
Arguments

x- a character or vector. FlyBase IDs or symbols.
n - The maximum number of genes to be opened. Default = 10.
```

Examples
```
j2g("e2f")
j2g(c("ovo", "otu"))
j2g("FBgn0000504")
j2g(c("FBgn0000504", "tra"))
```
