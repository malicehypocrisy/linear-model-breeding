# 一、分子血缘相关矩阵$\mathbf{A}$及其逆矩阵$\mathbf{A}^{-1}$的计算

## 1.$\mathbf{A}$的计算方法

### &#9312;通径系数

$$
a_{ii}=1+f_{i}\\
a_{ij}=\sum{(\frac{1}{2})^{n+n^{`}}(1+f_{A})}
$$

其中$f_{i}$为个体$i$的近交系数（为父母近交系数的一半），$n$和$n^`$是个体$i$和$j$经过不同通径到共同祖先$A$的时代数，$f_A$为$A$的近交系数。

### &#9313;循环公式

$$
a_{ii} = 1 + \frac{1}{2} a_{\text{父母}} \quad \text{（对角线元素）}
\\
a_{ij} = \frac{1}{2} \left( a_{i(j\text{父})} + a_{i(j\text{母})} \right) = a_{ji}
$$

> [!NOTE]
>
> **若个体$j$的父亲或母亲或双亲未知时，则$a_{ii}=1，a_{ij}=0$**

| 个体 |  父  |  母  |
| :--: | :--: | :--: |
|  1   |  -   |  -   |
|  2   |  -   |  -   |
|  3   |  1   |      |
|  4   |  1   |  2   |
|  5   |  3   |  4   |
|  6   |  1   |  4   |
|  7   |  5   |  6   |

列表注意：

&#9332;在个体一列中应包含所有在父和母列出现过的个体

&#9333;在个体一列保证后代不会出现在亲代前，一般按出生日期排序，先出生在前

&#9334;个体编号从1开始练习编号

## 2.$\mathbf{A}^{-1}$的计算方法

### &#9312;HENDERSON法

$$
\mathbf{A}^{-1}=(\mathbf{L}\mathbf{L}^{T})^{-1}=(\mathbf{TDD}\mathbf{T}^{T})^{-1}=\mathbf{T}^{T^{-1}}\mathbf{D}^{-2}\mathbf{T}^{-1}
$$

$\mathbf{T}^{-1}$（**下三角矩阵**）的对角线元素全为$1$，在其第$i$行个体的每个亲本对应的元素为-0.5，其余元素为零。

$\mathbf{D}^{-2}$是对角线矩阵$\delta_i$是其对角线上第$i$个元素
$$
\delta_i = \frac{1}{l_{ii}^2}\tag{1}
$$
**$l_{ii}$是$\mathbf{L}$阵的第$i$个对角线元素，$f$近交系数$f_i=a_{ii}-1$,$a_{ii}$是$\mathbf{A}$的对角线元素。**
$$
\\
l_{ii}=
\begin{cases}
[0.5 - 0.25(f_s + f_d)]^{0.5} & \text{当个体 i 的双亲 s 和 d 已知} \\
(0.75 - 0.25f_p)^{0.5} & \text{当个体 i 的一个亲体 p 已知} \\
1 & \text{当个体 i 双亲未知}\tag{2}
\end{cases}
$$

### &#9313;QUAAS法（简化$l_{ii}$的求解）

#### &#9372;**求解$l_{ii}$的步骤**：

&#9332;将所有个体按照计算$\mathbf{A}$阵时方式排列。

&#9333;建立两个阶数为n的零向量$\mathbf{v}$和$\mathbf{u}$,$\mathbf{v}$存放$l_{ii}$，$\mathbf{u}$用于存放$\mathbf{u}_i$,i=1,……，n,n为个体总数。

&#9334;对于i=1,……，n,计算

* $\mathbf{v}_i=l_{ii}(使用公式3)$

$$
u_i = \sum_{k=1}^i l_{ik}^2=a_{ii}
\\
l_{ii} =
\begin{cases}
\sqrt{1 - 0.25(u_s + u_d)} & \text{双亲 } s \text{ 和 } d \text{ 已知} \\
\sqrt{1 - 0.25u_p} & \text{单亲 } p \text{ 已知} \\
1 & \text{双亲未知}
\end{cases}\tag{3}
\\

\text{其中} (u_s, u_d, u_p) \text{为父母个体的累计值，反映遗传贡献}。
$$

* $\mathbf{u}_i=u_i+v_i^2$
* 对于i=1,……，n计算,

$$
l_{ki} =
\begin{cases}
\frac{l_{si} + l_{di}}{2} & \text{如个体 \(k\) 的双亲 \(s\) 和 \(d\) 排在 \(i\) 以后} \\
\frac{l_{pi}}{2} & \text{如个体 \(k\) 的一个亲体 \(p\) 排在 \(i\) 以后} \\
0 & \text{如个体 \(k\) 的双亲均排在 \(i\) 之前}
\end{cases}

\\

u_k = u_k + l_{ki}^2
$$

#### &#9373;求$\mathbf{A}^{-1}$的步骤

&#9332;将所有个按照$\mathbf{A}$求阵时的方式排列，

&#9333;将$\mathbf{A}^{-1}$中的所有元素置为零

&#9334;对于i=1,……，n,计算

* 按上述方法计算$l_{ii}$
* 如i的双亲s 和 d 已知，将下列数值加到$\mathbf{A}^{-1}$中：

$$
\begin{aligned}
\delta_i &\rightarrow a_{ii} \\
-\delta_i / 2 &\rightarrow a_{is}, a_{si}, a_{id}, a_{di} \\
\delta_i / 4 &\rightarrow a_{ss}, a_{dd}, a_{sd}, a_{ds}
\end{aligned}
$$

* 如i的一个亲本p已知，则
	$$
	\begin{aligned}
	\delta_i &\rightarrow a_{ii} \\
	-\delta_i / 2 &\rightarrow a_{ip}, a_{pi} \\
	\delta_i / 4 &\rightarrow a_{pp}
	\end{aligned}
	$$
	

* 如i的双亲均未知，则
	$$
	
	\delta_i \rightarrow a_{ii}
	$$
	

其中 $\delta_i$ 由式 (1) 计算，但当**群体为非近交群体时**，HENDERSON 证明：

$$

\delta_i =
\begin{cases}
2 & \text{当个体 } i \text{ 的双亲已知} \\
\frac{4}{3} & \text{当个体 } i \text{ 的一个亲体已知} \\
1 & \text{当个体 } i \text{ 双亲未知}
\end{cases} 
$$
