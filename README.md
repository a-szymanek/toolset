# toolset
Combines small tools that might be of use in the future

# Required setup

Before using tools user needs to specify the config file with connection string to the database. 

Create `development.ini` file in your home directory and add your connection string:

```bash
[Database]
sqlalchemy.url=postgres://username:password@hostname:port/databasename
```

Create `data` folder in the main directory
```bash
mkdir data
``` 


# List of tools

- germline_seq.py - downloading IMGT germline sequences with option of gapped/not gapped format

