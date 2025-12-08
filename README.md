Calculation and visualization of Hough modes
# Demo
https://nihsok.github.io/Hough/

Current setting
- Zonal wavenumber: s>0
- Latitudinal interval: 1 degree (fixed)
- Truncation number: 100 (under consideration)
- Treatment at singularity ($\tau$<12 hour): Interpolation

# English documentation
See [Wang et al. (2016)](https://doi.org/10.5194/gmd-9-1477-2016).

[Groves (1981)](https://doi.org/10.1016/0032-0633(81)90100-8) would also be helpful.

# Japanese documentation
## Wang et al. (2016)'s method
Hough関数$\Theta$は（規格化した）Legendre陪多項式$P_{r,s}$の重ね合わせとして得られる。

$\Theta=\sum_{r=s}^\infty a_r P_{r,s}(\mu)$

この形を潮汐方程式に代入することで、係数aの満たすべき式が得られる

$L_{r-2} a_{r-2}+(M_r-\lambda)a_r+L_r a_{r+2}=0$

これはaに対する連立方程式となっており、これを満たすaの組は行列の固有ベクトルとして求まる。

気温（およびジオポテンシャル、鉛直風）に関する$\Theta$が求まれば、東西風、南北風に関するHough関数もそこから計算できる。

この方法では周波数を固定して対応する等価深度を求めているが、等価深度を与えて周波数を計算する方法もある。