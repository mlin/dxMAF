# CrossSpeciesAlignments_Longest_Isoforms Developer Readme

<!--
TODO: Please edit this Readme.developer.md file to include information
for developers or advanced users, for example:

* Information about app internals and implementation details
* How to report bugs or contribute to development
-->

## Running this app with additional computational resources

This app has the following entry points:

* main
* process
* postprocess

When running this app, you can override the instance type to be used for each
entry point by providing the ``systemRequirements`` field to
```/applet-XXXX/run``` or ```/app-XXXX/run```, as follows:

    {
      systemRequirements: {
        "main": {"instanceType": "dx_m1.large"},
        "process": {"instanceType": "dx_m1.large"},
        "postprocess": {"instanceType": "dx_m1.large"}
      },
      [...]
    }

See <a
href="http://wiki.dnanexus.com/API-Specification-v1.0.0/IO-and-Run-Specifications#Run-Specification">Run
Specification</a> in the API documentation for more information about the
available instance types.
