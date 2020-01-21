# phylo

Computational Phylogenetics Tools in Rust

Note:

This code has a fair amount of functional use but is still undergoing cleanup in certain areas.  I needed different alignment handling
than rust-bio to handle certain PAML data, which is why this is currently implemented separately, though I will look again at that soon.  
I need to setup a fair amount of tests and complete missing chunks of the code before moving forward.

Also haven't fully setup Travis integration yet, or the external API.

Future:

#[![Build Status](https://travis-ci.org/den-sq/phylo.svg?branch=master)](https://travis-ci.org/den-sq/phylo)
#[![crates.io](https://img.shields.io/crates/v/phylo.svg)](https://crates.io/crates/phylo)


## Getting started

This project requires Rust to be installed. On OS X with Homebrew you can just run `brew install rust`.

Running it then should be as simple as:

```console
$ make
$ ./bin/phylo
```

### Testing

``make test``
