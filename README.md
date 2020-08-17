# Nanome - Molecular Dynamics [WIP - ALPHA Version]

Runs a simple molecular dynamics simulation on a complex of the workspace.

Needs a valid installation of OpenMM

![](molecular_dynamics.gif)

### Installation

```sh
$ pip install nanome-molecular-dynamics
```

### Usage

To start the plugin:

```sh
$ nanome-molecular-dynamics -a plugin_server_address
```

In Nanome:

- Activate Plugin, its window should open
- Select complex to use for simulation. Be careful: plugin will modify the complex to prepare it for simulation
- Click on Start

### Docker Usage

To run in a Docker container:

```sh
$ cd docker
$ ./build.sh
$ ./deploy.sh -a <plugin_server_address>
```

### License

MIT
