from aws_cdk import core

from clin_msi.stacks import ClinMsiStack


if __name__ == "__main__":
    app = core.App()
    ClinMsiStack(app, "clin-msi")
    app.synth()
