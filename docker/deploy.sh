if [ "$(docker ps -aq -f name=molecular-dynamics)" != "" ]; then
    # cleanup
    echo "removing exited container"
    docker rm -f molecular-dynamics
fi

docker run -d \
--name molecular-dynamics \
--restart unless-stopped \
-e ARGS="$*" \
molecular-dynamics
