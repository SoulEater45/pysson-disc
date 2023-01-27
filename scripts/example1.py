import plotly.graph_objects as go

from src.poisson_disc import poisson_disc_samples

width = 1920
height = 1080

samples = poisson_disc_samples(width, height, 10, 30)

print(samples.shape)
fig = go.Figure()
fig.add_trace(go.Scatter(x=samples[:, 0], y=samples[:, 1], mode='markers'))
fig.update_yaxes(
    scaleanchor = "x",
    scaleratio = 1,
)
fig.show()