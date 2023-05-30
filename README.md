# placenta
Workflow大部分在本地R studio中运行。  


SingleR注释在服务器上。  


# pesudotime given by monocle3效果欠佳  

√一种解释：the method（monocle） can only handle a single lineage？  

# 解决方案：

1、从fastq计算RNA velocity获得pesudo order；  

2、更换其他方法：

 Scanpy：Scanpy是一个强大的单细胞数据分析工具，它提供了一系列用于数据预处理、可视化、聚类、细胞发展轨迹分析等的函数和工具。它支持多种降维和聚类算法，并提供了多种细胞发展轨迹推断算法，如PAGA、DPT等。

Slingshot：Slingshot是一个用于单细胞RNA测序数据的细胞发展轨迹推断方法。它基于分歧时间点（branching events）来构建细胞发展轨迹，并可以用于多种细胞类型的数据分析。

Wanderlust：Wanderlust是另一个用于细胞发展轨迹推断的方法，它基于最小生成树（Minimum Spanning Tree）算法来构建细胞发展轨迹，并考虑了细胞类型的异质性。

TSCAN：TSCAN是一个用于细胞发展轨迹推断的方法，它利用时间信息和基因表达的变化模式来建模和推断细胞发展轨迹。

PAGA：PAGA（Partition-based Graph Abstraction）是一个基于图论的方法，用于可视化和推断细胞发展轨迹。它可以处理多种细胞类型的数据，并提供了图形化的展示和分析工具。
