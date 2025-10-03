#!/usr/bin/env python3
# plot_routes.py
# Uso: python3 plot_routes.py routes_data.txt
import sys
import os
import matplotlib.pyplot as plt
from matplotlib import cm

if len(sys.argv) < 2:
    print("Usage: python3 plot_routes.py <routes_data.txt>")
    sys.exit(1)

fname = sys.argv[1]
out_prefix = os.path.splitext(os.path.basename(fname))[0]

# parse file
with open(fname, "r") as f:
    first = f.readline().strip().split()
    # MUDANÇA: Agora lemos 4 valores do cabeçalho, incluindo a capacidade
    if len(first) < 4:
        raise RuntimeError("Formato inválido: header esperado: NPERIODS NCLIENTS NVEHICLES CAPACITY")
    Nperiods, nCust, nVeh, Capacity = map(int, first)
    
    # read depot line
    depot_line = f.readline().strip().split()
    if depot_line[0] != "DEPOT":
        raise RuntimeError("Expected 'DEPOT' line")
    depot_x = float(depot_line[1]); depot_y = float(depot_line[2])

    periods = []
    for _ in range(Nperiods):
        line = f.readline()
        while line and line.strip() == "":
            line = f.readline()
        if not line:
            break
        if not line.startswith("PERIOD"):
            raise RuntimeError("Expected 'PERIOD t' marker")
        
        header = f.readline().strip()
        if header != "CUSTOMERS":
            raise RuntimeError("Expected 'CUSTOMERS' after PERIOD")
        custs = []
        for _ in range(nCust):
            parts = f.readline().strip().split()
            if len(parts) < 7:
                raise RuntimeError("Linha de cliente com formato incorreto")
            custs.append({
                "id": int(parts[0]),
                "x": float(parts[1]), "y": float(parts[2]),
                "inv_before": float(parts[3]),
                "delivery": float(parts[4]),
                "demand": float(parts[5]),
                "inv_after": float(parts[6])
            })
        
        rtag = f.readline().strip()
        if rtag != "ROUTES":
            raise RuntimeError("Expected 'ROUTES' after CUSTOMERS block")
        routes = []
        while True:
            l = f.readline()
            if l is None:
                break
            l = l.strip()
            if l == "END_ROUTES":
                break
            if l == "":
                continue
            parts = l.split()
            route_nodes = [int(x) for x in parts]
            routes.append(route_nodes)
        periods.append({"customers": custs, "routes": routes})

# now plot per period
cmap = cm.get_cmap('tab10')
for t, per in enumerate(periods):
    fig, ax = plt.subplots(figsize=(9,6))
    custs = per["customers"]
    routes = per["routes"]

    # MUDANÇA: Cria um mapa de ID do cliente para a sua entrega para fácil acesso
    delivery_map = {c['id']: c['delivery'] for c in custs}

    # ... (código para plotar depot e clientes permanece o mesmo) ...
    ax.scatter([depot_x], [depot_y], marker="*", s=200, c="black", zorder=7)
    ax.annotate("DEPOT", xy=(depot_x, depot_y), xytext=(0,10), textcoords='offset points', ha='center', va='bottom', fontsize=9, bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.7), zorder=8)
    xs = [c["x"] for c in custs]
    ys = [c["y"] for c in custs]
    colors = ["#2ca02c" if c["inv_before"] - c["demand"] > 0 else ("#ffbf00" if abs(c["inv_before"] - c["demand"]) < 1e-9 else "#d62728") for c in custs]
    ax.scatter(xs, ys, s=100, c=colors, edgecolors='k', linewidths=0.6, zorder=6)
    for c in custs:
        ax.annotate(str(c["id"]), xy=(c["x"], c["y"]), xytext=(0,6), textcoords='offset points', ha='center', va='bottom', fontsize=8, bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.75), zorder=9, clip_on=False)

    # MUDANÇA: Lógica para desenhar as rotas foi atualizada
    for i, route in enumerate(routes):
        # Calcula a carga total da rota
        # O .get(node, 0) trata o caso do depot (node 0), que não tem entrega
        route_load = sum(delivery_map.get(node, 0) for node in route)
        
        # Define o estilo da linha com base na carga
        # Usamos uma pequena tolerância para casos de ponto flutuante
        if route_load < Capacity - 1e-9:
            linestyle = 'dotted'  # Pontilhado para rotas não cheias
        else:
            linestyle = 'solid'   # Sólido para rotas cheias

        # Converte nós para coordenadas
        coords = []
        for node in route:
            if node == 0:
                coords.append((depot_x, depot_y))
            else:
                # O ID do cliente é N, mas o índice na lista é N-1
                coords.append((custs[node - 1]["x"], custs[node - 1]["y"]))
        
        if len(coords) >= 2:
            xsr = [p[0] for p in coords]
            ysr = [p[1] for p in coords]
            color = cmap(i % 10)
            
            # Aplica o linestyle dinâmico
            ax.plot(xsr, ysr, linestyle=linestyle, linewidth=1.6, color=color, alpha=0.9, zorder=4)
            
            if len(xsr) > 2:
                ax.scatter(xsr[1:-1], ysr[1:-1], s=110, facecolors='none', edgecolors=[color], linewidths=1.6, zorder=7)

    ax.set_title(f"Período {t} — Rotas (Capacidade Máx. do Veículo: {Capacity})")
    ax.set_xlabel("X"); ax.set_ylabel("Y")
    ax.axis('equal')

    outname = f"{out_prefix}_period_{t}.png"
    plt.tight_layout()
    plt.savefig(outname, dpi=200)
    print(f"[plot_routes] saved {outname}")
    plt.close(fig)