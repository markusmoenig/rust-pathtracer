This is a port of the excellent [GLSL_Pathtracer](https://github.com/knightcrawler25/GLSL-PathTracer) to Rust utilizing an abstracted, trait based backend.

### Rust Features

* Multi-threading via [rayong]().
* An abstracted backend via a scene description trait. Implement hit tests, lights and materials for your procedural content.
* The pathtracer can be compiled to either ```f32``` or ```f64```. See the defines in the ```libs.rs```.

### How to use



### CPU vs GPU

