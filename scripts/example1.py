import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
from src.poisson_disc import poisson_disc_samples

width = 1920
height = 1080

bounds = np.array([width, height])
radii = np.array([10, 20, 30])
points, radii = poisson_disc_samples(bounds, radii, 30)

fig = go.Figure()
fig.update_yaxes(
    scaleanchor = "x",
    scaleratio = 1,
)
fig.add_trace(
    go.Scatter(
        x=points[:, 0],
        y=points[:, 1],
        mode='markers'
    )
)
for p, r in zip(points, radii):
    print(p[0], p[1], r)
    x0 = p[0]-r
    y0 = p[1]-r
    x1 = p[0]+r
    y1 = p[1]+r
    print(x0, y0, x1, y1)
    fig.add_shape(type="circle",
        xref="x", yref="y",
        x0=p[0]-r, y0=p[1]-r, x1=p[0]+r, y1=p[1]+r,
        line_color="green")

# fig.show()
pio.write_html(fig, file='out/example1.html')