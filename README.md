# minimum latency problem

## setup
```
git clone --recurse-submodules https://github.com/mororo18/mlp.git
```
or

```
git clone https://github.com/mororo18/mlp.git
cd mlp/
git submodule update --init
```
Then, you need to set up the [mlp-instances](https://github.com/mororo18/mlp-instances/blob/master/README.md) submodule repository.

## how to build (if necessary)

example:
```
cd java
./build.sh
```

## run

example:
```
cd mlp-instances/
./load att48.tsp
cd ../java
./run_java.sh
```
