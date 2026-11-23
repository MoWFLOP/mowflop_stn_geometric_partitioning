"""
RQ1 — disc_inv: invariância à discretização (Fig. 2).

Mantém o parque físico (mesmo sítio, mesmo τ) e refina a grade de posições
candidatas |P|; mede como o número de nós |N| da STN reage.
    Shannon (entropia, X = 70%) deve reagir: |N| cresce com |P|.
    Ocupação (grid, κ = 1,0) não deve: ℓ = κ·sqrt(A/τ) não contém |P|.

Usa as STNs agregadas (tags *norm, run = "0"), como nas demais figuras da RQ1.
Todas as instâncias de uma família devem ter o mesmo número de runs; o script
confere isso e avisa quando não for o caso. Para usar as STNs por execução
(mediana das runs), troque VARIANTE para "run".

disc_inv é a inclinação de log|N| contra log|P| (elasticidade):
    ≈ 0  → |N| não depende de |P| (invariante)
    ≈ 1  → |N| cresce na mesma proporção que |P|

Lê apenas as métricas já calculadas (pastas em ESQUEMAS):
    <campanha>/metrics/mowflop_<tag>norm_partition_metrics.csv     (VARIANTE = "agg")
    <campanha>/metrics/mowflop_<tag>runnorm_partition_metrics.csv  (VARIANTE = "run")
Quando as STNs das esparsas forem refeitas, basta trocar as pastas em ESQUEMAS.

Uso (a partir da pasta paper):
    python gerar_rq1_disc_inv.py
ou:
    python gerar_rq1_disc_inv.py caminho/para/paper

Saída: pasta rq1_disc_inv_saida/ com
    tabela_disc_inv.csv         |N| e |N₀| por instância, esquema e algoritmo
    disc_inv.csv                inclinação disc_inv por família, esquema e algoritmo
    grafico_disc_inv.png        Fig. 2: |N| x |P|, τ = 5 (ns506 a ns509)
    grafico_disc_inv_tau15.png  idem, τ = 15 (ns510 a ns513)
    artigo/                     versão para o artigo (IEEE, inglês), ver seção 4:
        disc_inv_tau5_moead.pdf   painel (a), MOEA/D
        disc_inv_tau5_nsga2.pdf   painel (b), NSGA-II
        disc_inv_tau5_legenda.pdf legenda única, para ir acima dos dois painéis
        (e os mesmos nomes em .png, 600 dpi, para conferência)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
OUT = ROOT / "rq1_disc_inv_saida"
OUT.mkdir(exist_ok=True)

P = 100          # configuração p100i50
VARIANTE = "agg"  # "agg": STN agregada | "run": mediana das STNs por execução

# Famílias com o mesmo parque físico e o mesmo τ; só |P| varia.
FAMILIAS = {
    5: ["ns506", "ns507", "ns508", "ns509"],
    15: ["ns510", "ns511", "ns512", "ns513"],
}

# (pasta da campanha, tag sem o sufixo norm/runnorm, rótulo, cor, marcador)
ESQUEMAS = [
    ("../STNs-MOCO-MoWFLOP/sharing/campanhas/full_ext_entropy_sparse_2026-10-05", "x70", "Shannon, X = 70%", "tab:blue", "o"),
    ("../STNs-MOCO-MoWFLOP/sharing/campanhas/full_ext_grid_sparse_2026-10-05", "g1.0", "Ocupação (grid), κ = 1,0", "tab:orange", "s"),
]
SUFIXO = {"agg": "norm", "run": "runnorm"}[VARIANTE]
ROT_N = {"agg": "|N| (STN agregada)", "run": "|N| (mediana das STNs por execução)"}[VARIANTE]
ALGOS = [("MOEAD", "MOEA/D"), ("NSGA2", "NSGA-II")]


def ler(campanha, tag):
    f = ROOT / campanha / "metrics" / f"mowflop_{tag}{SUFIXO}_partition_metrics.csv"
    d = pd.read_csv(f, dtype={"run": str})
    return d[d["p"] == P]


# ---------------------------------------------------------------------------
# 0) Conferência: mesmo número de runs em todas as instâncias de cada família
# ---------------------------------------------------------------------------
for campanha, tag, rotulo, _, _ in ESQUEMAS:
    f = ROOT / campanha / "metrics" / f"mowflop_{tag}runnorm_partition_metrics.csv"
    if not f.exists():
        print(f"[aviso] {f.name} não encontrado; número de runs não conferido ({rotulo})")
        continue
    r = pd.read_csv(f, dtype={"run": str})
    r = r[r["p"] == P].groupby(["instance", "algo"]).size()
    for tau, insts in FAMILIAS.items():
        n = r[r.index.get_level_values("instance").isin(insts)]
        if n.nunique() > 1:
            print(f"[aviso] {rotulo}, τ = {tau}: número de runs difere entre instâncias:")
            print(n.unstack("algo").to_string())
        elif len(n):
            print(f"{rotulo}, τ = {tau}: {n.iloc[0]} runs em todas as instâncias")


# ---------------------------------------------------------------------------
# 1) |N| (e |N₀|, sem partição) por instância e algoritmo
# ---------------------------------------------------------------------------
esparsas = sum(FAMILIAS.values(), [])
partes = []
for campanha, tag, rotulo, _, _ in ESQUEMAS:
    d = ler(campanha, tag)
    d = d[d["instance"].isin(esparsas)]
    if VARIANTE == "agg":
        t = d[["instance", "tau", "n_positions", "algo", "nodes", "nodes_raw"]].copy()
    else:
        t = (d.groupby(["instance", "tau", "n_positions", "algo"])
              .agg(runs=("run", "size"), nodes=("nodes", "median"), nodes_raw=("nodes_raw", "median"))
              .reset_index())
    t["esquema"] = rotulo
    partes.append(t)
tab = pd.concat(partes).sort_values(["tau", "esquema", "algo", "n_positions"])
tab.to_csv(OUT / "tabela_disc_inv.csv", index=False)
print(f"{ROT_N}, p = {P}:")
print(tab.to_string(index=False))

# ---------------------------------------------------------------------------
# 2) disc_inv = inclinação de log10|N| x log10|P|
# ---------------------------------------------------------------------------
def inclinacao(x, y):
    return np.polyfit(np.log10(x), np.log10(y), 1)[0]

linhas = []
for tau in FAMILIAS:
    for (esq, algo), g in tab[tab["tau"] == tau].groupby(["esquema", "algo"]):
        linhas.append({"tau": tau, "esquema": esq, "algo": algo,
                       "instancias": ", ".join(g["instance"]),
                       "disc_inv": round(inclinacao(g["n_positions"], g["nodes"]), 3)})
    # referência: STN sem partição (|N₀| é o mesmo nos dois esquemas)
    for algo, g in tab[(tab["tau"] == tau) & (tab["esquema"] == ESQUEMAS[0][2])].groupby("algo"):
        linhas.append({"tau": tau, "esquema": "sem partição |N₀|", "algo": algo,
                       "instancias": ", ".join(g["instance"]),
                       "disc_inv": round(inclinacao(g["n_positions"], g["nodes_raw"]), 3)})
di = pd.DataFrame(linhas)
di.to_csv(OUT / "disc_inv.csv", index=False)
print("\ndisc_inv (inclinação de log|N| x log|P|):")
print(di.to_string(index=False))

# ---------------------------------------------------------------------------
# 3) Figura: |N| x |P|, uma coluna por algoritmo
# ---------------------------------------------------------------------------
for tau, insts in FAMILIAS.items():
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, (algo, rot_alg) in zip(axs, ALGOS):
        t = tab[(tab["tau"] == tau) & (tab["algo"] == algo)]
        ref = t[t["esquema"] == ESQUEMAS[0][2]]
        ax.plot(ref["n_positions"], ref["nodes_raw"], "^--", color="0.6", label="sem partição |N₀|")
        for _, _, rotulo, cor, mk in ESQUEMAS:
            s = t[t["esquema"] == rotulo]
            ax.plot(s["n_positions"], s["nodes"], marker=mk, color=cor, label=rotulo)
            for _, r in s.iterrows():
                # rótulo acima da curva de cima e abaixo da de baixo, para não sobrepor
                outro = t[(t["instance"] == r["instance"]) & (t["esquema"] != rotulo)]["nodes"]
                dy = 4 if outro.empty or r["nodes"] >= outro.iloc[0] else -12
                ax.annotate(r["instance"], (r["n_positions"], r["nodes"]), fontsize=8,
                            xytext=(4, dy), textcoords="offset points", color=cor)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("posições candidatas |P|")
        ax.set_title(f"{rot_alg}, τ = {tau}, p = {P}")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
    axs[0].set_ylabel(ROT_N)
    fig.tight_layout()
    nome = "grafico_disc_inv.png" if tau == 5 else f"grafico_disc_inv_tau{tau}.png"
    fig.savefig(OUT / nome, dpi=200)
    plt.close(fig)

print(f"\nArquivos gravados em {OUT}")

# ---------------------------------------------------------------------------
# 4) Figuras para o artigo (IEEE, duas colunas)
# ---------------------------------------------------------------------------
# Regras dos revisores:
#   - fonte da figura do mesmo tamanho que a do texto: a figura é gerada já na
#     largura impressa (sem redução no LaTeX) e com fonte de 8 pt;
#   - sem título no gráfico: o algoritmo vai no subfloat (a)/(b) do LaTeX;
#   - legenda sempre na mesma ordem; aqui ela é única, num arquivo à parte,
#     posta acima dos dois painéis.
# LaTeX sugerido (sem nenhum "width" diferente destes, para não reescalar):
#   \includegraphics{disc_inv_tau5_legenda.pdf}\\
#   \subfloat[MOEA/D]{\includegraphics{disc_inv_tau5_moead.pdf}}\hfill
#   \subfloat[NSGA-II]{\includegraphics{disc_inv_tau5_nsga2.pdf}}
ART = OUT / "artigo"
ART.mkdir(exist_ok=True)

COLUNA_IN = 3.5                 # largura da coluna IEEE (\columnwidth = 252 pt)
PAINEL_IN = 0.49 * COLUNA_IN    # cada subfloat ocupa 0.49\columnwidth
ALTURA_IN = 1.45
FONTE_PT = 8                    # caption IEEE = 8 pt, texto = 10 pt
FAMILIAS_ARTIGO = [5]           # só τ = 5 entra no artigo

# Rótulos em inglês, na ordem em que aparecem na legenda (e no texto).
# Os parâmetros (X = 70%, κ = 1) ficam no texto do artigo, não na figura.
ROT_ARTIGO = {
    "ref": "Unpartitioned",
    ESQUEMAS[0][2]: "Shannon",
    ESQUEMAS[1][2]: "Occupancy grid",
}

estilo = {
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Liberation Serif", "Nimbus Roman",
                   "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": FONTE_PT,
    "axes.labelsize": FONTE_PT,
    "xtick.labelsize": FONTE_PT,
    "ytick.labelsize": FONTE_PT,
    "legend.fontsize": FONTE_PT,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "lines.linewidth": 1.0,
    "lines.markersize": 3.5,
    "pdf.fonttype": 42,          # fontes embutidas (exigência do PDF eXpress)
    "ps.fonttype": 42,
}


def series(t):
    """Séries na ordem fixa da legenda: sem partição, Shannon, grid."""
    ref = t[t["esquema"] == ESQUEMAS[0][2]]
    yield ("ref", ref["n_positions"], ref["nodes_raw"],
           dict(color="0.55", marker="^", linestyle="--"))
    for _, _, rotulo, cor, mk in ESQUEMAS:
        s = t[t["esquema"] == rotulo]
        yield (rotulo, s["n_positions"], s["nodes"],
               dict(color=cor, marker=mk, linestyle="-"))


def salvar(fig, nome):
    for ext, kw in ((".pdf", {}), (".png", {"dpi": 600})):
        fig.savefig(ART / f"{nome}{ext}", **kw)
    plt.close(fig)


with plt.rc_context(estilo):
    for tau in FAMILIAS_ARTIGO:
        t_tau = tab[tab["tau"] == tau]
        # mesmos limites nos dois painéis, para que sejam comparáveis
        ys = pd.concat([t_tau["nodes"], t_tau["nodes_raw"]])
        ylim = (ys.min() / 1.6, ys.max() * 1.6)
        xs = t_tau["n_positions"]
        xlim = (xs.min() / 1.6, xs.max() * 1.6)

        for algo, _ in ALGOS:
            t = t_tau[t_tau["algo"] == algo]
            fig = plt.figure(figsize=(PAINEL_IN, ALTURA_IN))
            # margens fixas: os dois painéis ficam com eixos do mesmo tamanho
            ax = fig.add_axes([0.25, 0.25, 0.71, 0.71])
            for chave, x, y, kw in series(t):
                ax.plot(x, y, **kw)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel(r"$|\mathcal{P}|$", labelpad=1)
            ax.set_ylabel(r"$|\mathcal{N}|$", labelpad=1)
            ax.tick_params(pad=1.5)
            ax.grid(alpha=0.3, which="major", linewidth=0.4)
            salvar(fig, f"disc_inv_tau{tau}_{algo.lower()}")

        # legenda única, uma linha, na largura da coluna
        fig = plt.figure(figsize=(COLUNA_IN, 0.22))
        handles = [plt.Line2D([], [], **kw) for _, _, _, kw in series(t_tau)]
        labels = [ROT_ARTIGO[k] for k, _, _, _ in series(t_tau)]
        fig.legend(handles, labels, loc="center", ncol=len(handles), frameon=False,
                   handlelength=2.2, columnspacing=1.2, handletextpad=0.5,
                   borderaxespad=0, borderpad=0)
        salvar(fig, f"disc_inv_tau{tau}_legenda")

print(f"Figuras do artigo em {ART}")

