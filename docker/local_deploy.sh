if [ "$(docker ps -aq -f name=molecular-dynamics)" != "" ]; then
        # cleanup
        echo "removing exited container"
        docker rm -f molecular-dynamics
fi

if [ "$1" != "" ]; then
    echo "Using specified plugin server: $1"
    docker run -d \
    -p 8888:8888 \
    -e PLUGIN_SERVER=$1 \
    --name molecular-dynamics molecular-dynamics
else
    echo "Using default plugin server: plugins.nanome.ai"
    docker run -d \
    -p 8888:8888 \
    --name molecular-dynamics molecular-dynamics
fi